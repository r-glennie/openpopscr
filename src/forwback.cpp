// Copyright (c) 2019 Richard Glennie, University of St Andrews
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
// moveds.cpp: C++ functions to compute conditional probabilities 
//
////////////////////////////////////////////////////////////////////////////////
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>

using namespace RcppParallel; 

struct AlphaCalculator : public Worker {
  
  // input 
  const int n;
  const int J; 
  const int M; 
  const arma::mat pr0; 
  const Rcpp::List pr_capture;
  const Rcpp::List tpms;
  const int num_states;
  const arma::vec entry;
  
  // transform 
  std::vector<arma::mat> tpm; 
  std::vector<arma::cube> pr_cap; 

  // output 
  arma::field<arma::cube>& lalpha; 
  
  // initialiser
  AlphaCalculator(const int n, const int J, const int M, 
                const arma::mat pr0, 
                const Rcpp::List pr_capture, 
                const Rcpp::List tpms,
                const int num_states,
                const arma::vec entry,
                arma::field<arma::cube>& lalpha) : n(n), J(J), M(M), pr0(pr0), pr_capture(pr_capture), tpms(tpms), num_states(num_states), entry(entry), lalpha(lalpha) {
    if (num_states > 1) {
      tpm.resize(J); 
      for (int j = 0; j < J - 1; ++j) tpm[j] = Rcpp::as<arma::mat>(tpms[j]); 
    }
    pr_cap.resize(n);
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector pr_capvec(pr_capture[i]);
      arma::cube pr_icap(pr_capvec.begin(), M, num_states, J, false);
      pr_cap[i] = pr_icap;
    }
  }
 
  void operator()(std::size_t begin, std::size_t end) { 
    
    for (int i = begin; i < end; ++i) {
      double llk = 0; 
      double sum_pr; 
      arma::mat pr = pr0; 
      arma::cube prcap; 
      for (int j = entry(i); j < J - 1; ++j) {
        pr %= pr_cap[i].slice(j);
        if (num_states > 1) {
          pr *= tpm[j]; 
        }
        sum_pr = accu(pr); 
        llk += log(sum_pr); 
        pr /= sum_pr; 
        lalpha(i).slice(j) = log(pr) + llk; 
      }
      pr %= pr_cap[i].slice(J - 1);
      sum_pr = accu(pr); 
      llk += log(sum_pr); 
      pr /= sum_pr; 
      lalpha(i).slice(J - 1) = log(pr) + llk; 
    }
  }
};


//' Computes forward probabilities
//'
//' @param n number of individuals 
//' @param J total number of occasions 
//' @param M total number of mesh points
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_states 2 = CJS model, 3 = JS model 
//' @param entry vector of entry occasions per individual 
//'
//' @return log-likelihood value 
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_alpha(const int n, const int J, const int M, 
                  const arma::mat pr0, 
                  const Rcpp::List pr_capture, 
                  const Rcpp::List tpms,
                  const int num_states,
									const arma::vec entry) {
  
  arma::field<arma::cube> lalpha(n);
  for (int i = 0; i < n; ++i) lalpha(i) = arma::zeros<arma::cube>(M, num_states, J); 
  AlphaCalculator alpha_calculator(n, J, M, pr0, pr_capture, tpms, num_states, entry, lalpha); 
  parallelFor(0, n, alpha_calculator, 1); 
  return(lalpha); 
}

struct BetaCalculator : public Worker {
  
  // input 
  const int n;
  const int J; 
  const int M; 
  const arma::mat pr0; 
  const Rcpp::List pr_capture;
  const Rcpp::List tpms;
  const int num_states;
  const arma::vec entry;
  
  // transform 
  std::vector<arma::mat> tpm; 
  std::vector<arma::cube> pr_cap; 
  
  // output 
  arma::field<arma::cube>& lbeta; 
  
  // initialiser
  BetaCalculator(const int n, const int J, const int M, 
                  const arma::mat pr0, 
                  const Rcpp::List pr_capture, 
                  const Rcpp::List tpms,
                  const int num_states,
                  const arma::vec entry,
                  arma::field<arma::cube>& lbeta) : n(n), J(J), M(M), pr0(pr0), pr_capture(pr_capture), tpms(tpms), num_states(num_states), entry(entry), lbeta(lbeta) {
    if (num_states > 1) {
      tpm.resize(J); 
      for (int j = 0; j < J - 1; ++j) tpm[j] = Rcpp::as<arma::mat>(tpms[j]); 
    }
    pr_cap.resize(n);
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector pr_capvec(pr_capture[i]);
      arma::cube pr_icap(pr_capvec.begin(), M, num_states, J, false);
      pr_cap[i] = pr_icap;
    }
  }
  
  void operator()(std::size_t begin, std::size_t end) { 
    
    for (int i = begin; i < end; ++i) {
      arma::mat pr(M, num_states, arma::fill::ones);
      pr /= (1.0 * M * num_states); 
      double llk = log((1.0 * M * num_states)); 
      double sum_pr; 
      arma::cube prcap; 
      lbeta(i).slice(J - 1).zeros(); 
      for (int j =  J - 2; j > entry(i) - 1; --j) {
        pr %= pr_cap[i].slice(j + 1);
        if (num_states > 1) {
          pr = tpm[j] * pr.t(); 
        }
        lbeta(i).slice(j) = log(pr) + llk; 
        sum_pr = accu(pr); 
        llk += log(sum_pr); 
        pr /= sum_pr; 
      }
    }
  }
};


//' Computes backward probabilities
//'
//' @param n number of individuals 
//' @param J total number of occasions 
//' @param M total number of mesh points
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_states 2 = CJS model, 3 = JS model 
//' @param entry vector of entry occasions per individual 
//'
//' @return log-likelihood value 
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_beta(const int n, const int J, const int M, 
                                     const arma::mat pr0, 
                                     const Rcpp::List pr_capture, 
                                     const Rcpp::List tpms,
                                     const int num_states,
                                     const arma::vec entry) {
  
  arma::field<arma::cube> lbeta(n);
  for (int i = 0; i < n; ++i) lbeta(i) = arma::zeros<arma::cube>(M, num_states, J); 
  BetaCalculator beta_calculator(n, J, M, pr0, pr_capture, tpms, num_states, entry, lbeta); 
  parallelFor(0, n, beta_calculator, 1); 
  return(lbeta); 
}
