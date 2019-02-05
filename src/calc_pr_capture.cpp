// Copyright (c) 2017 Richard Glennie, University of St Andrews
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files, to deal in the software
// without restriction, including without limitation the right to use, copy,
// modify, publish, distribute, sublicense, and/or sell copies of the software,
// and to permit persons to whom the software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS  OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE THE USE OR OTHER DEALINGS IN
// THE SOFTWARE
//
////////////////////////////////////////////////////////////////////////////////
// openpopscr project: open population spatial capture-recapture modelling
//
// moveds.cpp: C++ functions to compute capture history probabilities 
//
////////////////////////////////////////////////////////////////////////////////
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel; 

struct PrCaptureCalculator : public Worker {
  
  // input 
  const int J; 
  const int K; 
  const int M;
  const int alive_col; 
  const arma::cube capthist;
  const arma::cube enc0;
  const arma::mat usage;
  const int num_states;
  const int detector_type;
  const int n_prim; 
  const arma::vec S; 
  
  // output 
  arma::field<arma::cube>& probfield;
  
  // initialise
  PrCaptureCalculator(const int J, 
                        const int K, 
                        const int M, 
                        const int alive_col, 
                        const arma::cube capthist, 
                        const arma::cube enc0, 
                        const arma::mat usage, 
                        const int num_states,
                        const int detector_type,
                        const int n_prim, 
                        const arma::vec S, 
                        arma::field<arma::cube>& probfield) : J(J), K(K), M(M), 
                        alive_col(alive_col), 
                        capthist(capthist), 
                        enc0(enc0), 
                        usage(usage), 
                        num_states(num_states), 
                        detector_type(detector_type),
                        n_prim(n_prim), 
                        S(S), 
                        probfield(probfield) {} 
  
  void operator()(std::size_t begin, std::size_t end) { 
    arma::cube iprob = arma::zeros<arma::cube>(M, num_states, n_prim);
    arma::vec enc; 
    arma::vec penc; 
    arma::vec savedenc;
    arma::mat encslice; 
    double num; 
    arma::vec probslice = arma::zeros<arma::vec>(M);
    double sumcap; 
    for (int i = begin; i < end; ++i) {
      int j = -1; 
      for (int prim = 0; prim < n_prim; ++prim) { 
        probslice.zeros(); 
        bool unseen = true;
        for (int s = 0; s < S(prim); ++s) {
          ++j; 
          encslice = enc0.slice(j);
          penc = arma::zeros<arma::vec>(M); 
          savedenc = arma::zeros<arma::vec>(M); 
          sumcap = 0; 
          for (int k = 0; k < K; ++k) {
            if (usage(k, j) < 1e-16) continue; 
            enc = encslice.col(k) * usage(k, j); 
            if (detector_type == 1 | detector_type == 4) {
              // avoid zeros 
              enc += 1e-16;
              probslice += capthist(i, j, k) * log(enc) - enc;// - lgamma(capthist(i, j, k) + 1); 
            } else if (detector_type == 2) {
              penc = 1.0 - exp(-enc); 
              // avoid zeros 
              penc += 1e-16; 
              probslice += capthist(i, j, k) * log(penc) - (1.0 - capthist(i, j, k)) * enc; 
            } else if (detector_type == 3) {
              penc += enc; 
              sumcap += capthist(i, j, k); 
              if(capthist(i, j, k) > 0.5) savedenc += capthist(i, j, k) * log(enc); 
            }
            if (capthist(i, j, k) > 1e-16) unseen = false; 
          }
          if (detector_type == 3) {
            // avoid zeros
            penc += 1e-16; 
            if (!unseen) probslice += savedenc - sumcap * log(penc); 
            probslice += -(1.0 - sumcap) * penc + sumcap * log(1.0 - exp(-penc)); 
          }
        }
        if (unseen) {
          if (num_states == 2) iprob.slice(prim).col(1 - alive_col).ones(); 
          if (num_states == 3) {
            iprob.slice(prim).col(0).ones(); 
            iprob.slice(prim).col(2).ones();  
          }
        }
        iprob.slice(prim).col(alive_col) = exp(probslice); 
      }
     probfield(i) = iprob;
    }
  }
};  
  
//' Computes probability of each capture record
//'
//' @param n number of individuals 
//' @param J total number of occasions 
//' @param K total number of traps ever used  
//' @param M total number of mesh points
//' @param capvec a pointer to the capthist array 
//' @param enc_rate a pointer to the encounter rate array, see calc_pr_capture() in JsModel
//' @param usage matrix with J x K where (j,k) entry is usage of trap k in occasion j
//' @param num_cores number of processor cores to use in parallelisation 
//' @param num_states: 1 = SCR model, 2 = CJS model, 3 = JS model 
//' @param detector_type 1 = count, 2 = proximity/binary, 3 = multi-catch, 4 = transect 
//'
//' @return  Array with (i,j,m) entry the probability of capture record for individual i in occasion j given activity centre at mesh point m  
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_pr_capture(const int n, const int J, const int K, const int M, 
                             Rcpp::NumericVector& capvec, 
                             Rcpp::NumericVector& enc_rate,
                             const arma::mat usage, 
                             const int num_cores, 
                             const int num_states,
                             const int detector_type, 
                             const int n_prim, 
                             const arma::vec S) {
  
  const arma::cube capthist(capvec.begin(), n, J, K, false);
  const arma::cube enc0(enc_rate.begin(), M, K, J, false);
  int alive_col = 1; 
  if (num_states < 3) alive_col = 0; 
  arma::field<arma::cube> probfield(n);
  PrCaptureCalculator pr_capture_calc(J, K, M, alive_col, capthist, enc0, usage, num_states, detector_type, n_prim, S, probfield); 
  parallelFor(0, n, pr_capture_calc); 
  return(probfield);
}
