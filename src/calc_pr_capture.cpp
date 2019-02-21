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
  const arma::cube& capthist;
  const arma::cube& enc0;
  const arma::mat& usage;
  const int num_states;
  const int detector_type;
  const int n_prim; 
  const arma::vec& S; 
  const arma::vec& entry; 
  arma::mat total_enc;
  arma::mat log_total_enc; 
  arma::mat log_total_penc; 
  arma::cube logenc0; 
  arma::cube log_penc; 

  // output 
  arma::field<arma::cube>& probfield;
  
  // initialise
  PrCaptureCalculator(const int J, 
                        const int K, 
                        const int M, 
                        const int alive_col, 
                        const arma::cube& capthist, 
                        const arma::cube& enc0, 
                        const arma::mat& usage, 
                        const int num_states,
                        const int detector_type,
                        const int n_prim, 
                        const arma::vec& S, 
                        const arma::vec& entry, 
                        arma::field<arma::cube>& probfield) : J(J), K(K), M(M), 
                        alive_col(alive_col), 
                        capthist(capthist), 
                        enc0(enc0), 
                        usage(usage), 
                        num_states(num_states), 
                        detector_type(detector_type),
                        n_prim(n_prim), 
                        S(S), 
                        entry(entry), 
                        probfield(probfield) {
    
    if (detector_type == 3) { 
      total_enc = arma::zeros<arma::mat>(M, J); 
      log_total_enc = arma::zeros<arma::mat>(M, J); 
      log_total_penc = arma::zeros<arma::mat>(M, J); 
      int j = -1; 
      for (int prim = 0; prim < n_prim; ++prim) {
        for (int s = 0; s < S(prim); ++s) {
          ++j; 
          total_enc.col(j) = enc0.slice(j) * usage.col(j); 
        }
      }
      log_total_enc = log(total_enc); 
      log_total_penc = 1.0 - exp(-total_enc); 
      log_total_penc = log(log_total_penc + 1e-16); 
    } else if (detector_type == 2) {
      log_penc = arma::zeros<arma::cube>(M, K, J); 
      int j = -1; 
      for (int prim = 0; prim < n_prim; ++prim) {
        for (int s = 0; s < S(prim); ++s) {
          ++j; 
          for (int k = 0; k < K; ++k) { 
            log_penc.slice(j).col(k) = 1.0 - exp(-enc0.slice(j).col(k) * usage(k, j)); 
            log_penc.slice(j).col(k) = log(log_penc.slice(j).col(k) + 1e-16); 
          }
        }
      }
    }
    logenc0 = log(enc0); 
  } 
  
  void operator()(std::size_t begin, std::size_t end) { 
    for (int i = begin; i < end; ++i) {
      arma::vec savedenc(M);
      double sumcap;
      int j = -1; 
      for (int prim = 0; prim < n_prim; ++prim) { 
       bool unseen = true;
        if (entry(i) - 1 < prim) {
         for (int s = 0; s < S(prim); ++s) {
           ++j;
           if (detector_type == 3) savedenc.zeros();
           sumcap = 0;
            for (int k = 0; k < K; ++k) {
              if (usage(k, j) < 1e-16) continue;
              if (capthist(i, j, k) < 0.5 & detector_type == 3) continue;
              if (detector_type == 1 | detector_type == 4) {
                probfield(i).slice(prim).col(alive_col) += capthist(i, j, k) * logenc0.slice(j).col(k) - enc0.slice(j).col(k);
              } else if (detector_type == 2) {
                probfield(i).slice(prim).col(alive_col) += capthist(i, j, k) * log_penc.slice(j).col(k) - (1.0 - capthist(i, j, k)) * enc0.slice(j).col(k);
              } else if (detector_type == 3) {
                sumcap += capthist(i, j, k);
                savedenc += capthist(i, j, k) * logenc0.slice(j).col(k);
              }
              if (capthist(i, j, k) > 1e-16) unseen = false;
            }
            if (detector_type == 3) {
              if (!unseen) probfield(i).slice(prim).col(alive_col) += savedenc - sumcap * log_total_enc.col(j); 
              probfield(i).slice(prim).col(alive_col) += -(1.0 - sumcap) * total_enc.col(j) + sumcap * log_total_penc.col(j); 
            }
          }
        }
        if (unseen) {
          if (num_states == 2) probfield(i).slice(prim).col(1 - alive_col).ones(); 
          if (num_states == 3) {
            probfield(i).slice(prim).col(0).ones();
            probfield(i).slice(prim).col(2).ones();
          }
        }
        probfield(i).slice(prim).col(alive_col) = exp(probfield(i).slice(prim).col(alive_col)); 
      }
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
//' @param num_states: 1 = SCR model, 2 = CJS model, 3 = JS model 
//' @param detector_type 1 = count, 2 = proximity/binary, 3 = multi-catch, 4 = transect 
//'
//' @return  Array with (i,j,m) entry the probability of capture record for individual i in occasion j given activity centre at mesh point m  
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_pr_capture(const int n, const int J, const int K, const int M, 
                             const arma::cube& capthist, 
                             const arma::cube& enc0,
                             const arma::mat usage, 
                             const int num_states,
                             const int detector_type, 
                             const int n_prim, 
                             const arma::vec S, 
                             const arma::vec entry) {
  int alive_col = 1; 
  if (num_states < 3) alive_col = 0; 
  arma::field<arma::cube> probfield(n);
  for (int i = 0; i < n; ++i) probfield(i) = arma::zeros<arma::cube>(M, num_states, n_prim); 
  PrCaptureCalculator pr_capture_calc(J, K, M, alive_col, capthist, enc0, usage, num_states, detector_type, n_prim, S, entry, probfield); 
  parallelFor(0, n, pr_capture_calc); 
  return(probfield);
}
