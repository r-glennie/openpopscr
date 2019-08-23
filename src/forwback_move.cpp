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
#include "forward_declare.h"

using namespace RcppParallel; 

struct AlphaMoveCalculator : public Worker {
  
  // input 
  const int n; 
  const int J;
  const arma::mat pr0; 
  const Rcpp::List pr_capture; 
  const Rcpp::List tpms;
  const arma::vec num_cells; 
  const arma::vec inside; 
  const double dx; 
  const arma::vec dt; 
  const arma::mat sd; 
  const int num_states;
  const int minstate; 
  const int maxstate; 
  const arma::vec entry; 
  
  // transform 
  std::vector<arma::mat> tpm; 
  std::vector<arma::cube> pr_cap; 
  std::vector<arma::sp_mat> trm; 
  
  // output 
  arma::field<arma::cube>& lalpha; 
  
  // initialiser
  AlphaMoveCalculator(const int n, const int J, 
                      const arma::mat pr0, 
                      const Rcpp::List pr_capture, 
                      const Rcpp::List tpms,
                      const arma::vec num_cells, 
                      const arma::vec inside, 
                      const double dx, 
                      const arma::vec dt, 
                      const arma::mat sd,
                      const int num_states,
                      const int minstate, 
                      const int maxstate, 
                      const arma::vec entry,
                      arma::field<arma::cube>& lalpha) : n(n), J(J), pr0(pr0), pr_capture(pr_capture), tpms(tpms), num_cells(num_cells), inside(inside), dx(dx), dt(dt), sd(sd), num_states(num_states), minstate(minstate), maxstate(maxstate), entry(entry), lalpha(lalpha) {
    if (num_states > 1) {
      tpm.resize(J); 
      for (int j = 0; j < J - 1; ++j) tpm[j] = Rcpp::as<arma::mat>(tpms[j]); 
    }
    trm.resize(J * num_states); 
    for (int j = 0; j < J - 1; ++j) {
      for (int g = minstate; g < minstate + num_states; ++g) {
        if (sd(j, g - minstate) < 0) continue; 
        trm[g - minstate + j * num_states] = CalcTrm(num_cells, sd(j, g - minstate), dx, inside); 
      }
    }
    pr_cap.resize(n);
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector pr_capvec(pr_capture[i]);
      arma::cube pr_icap(pr_capvec.begin(), num_cells(0), num_states, J, false);
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
        for (int g = minstate; g < minstate + num_states; ++g) {
          if (sd(j, g - minstate) < 0) continue; 
          try {
            pr.col(g) = ExpG(pr.col(g), trm[g - minstate + j * num_states], dt(j));
          } catch(...) {
            break;
          }
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
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_cells number of cells in x,y,total 
//' @param inside 0 if meshpt outside survey region, 1 otherwise 
//' @param dx mesh spacing 
//' @param dt time between occasions 
//' @param sd movement parameter for each occasion 
//' @param num_states 2 = CJS model, 3 = JS model 
//' @param entry time each individual entered survey 
//'
//' @return forwards probabiltiies individual x occasion x mesh x state
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_movealpha(const int n, const int J,
                                         const arma::mat pr0, 
                                         const Rcpp::List pr_capture, 
                                         const Rcpp::List tpms,
                                         const arma::vec num_cells, 
                                         const arma::vec inside, 
                                         const double dx, 
                                         const arma::vec dt, 
                                         const arma::mat sd, 
                                         const int num_states,
                                         const int minstate, 
                                         const int maxstate, 
                                         const arma::vec entry) {
  
  arma::field<arma::cube> lalpha(n);
  for (int i = 0; i < n; ++i) lalpha(i) = arma::zeros<arma::cube>(num_cells(0), num_states, J); 
  AlphaMoveCalculator alpha_calculator(n, J, pr0, pr_capture, tpms, num_cells, inside, dx, dt, sd, num_states, minstate, maxstate, entry, lalpha); 
  parallelFor(0, n, alpha_calculator, 1); 
  return(lalpha); 
}

struct BetaMoveCalculator : public Worker {
  
  // input 
  const int n; 
  const int J;
  const arma::mat pr0; 
  const Rcpp::List pr_capture; 
  const Rcpp::List tpms;
  const arma::vec num_cells; 
  const arma::vec inside; 
  const double dx; 
  const arma::vec dt; 
  const arma::mat sd; 
  const int num_states;
  const int minstate; 
  const int maxstate; 
  const arma::vec entry; 
  
  // transform 
  std::vector<arma::mat> tpm; 
  std::vector<arma::cube> pr_cap; 
  std::vector<arma::sp_mat> trm; 
  
  // output 
  arma::field<arma::cube>& lbeta; 
  
  // initialiser
  BetaMoveCalculator(const int n, const int J, 
                     const arma::mat pr0, 
                     const Rcpp::List pr_capture, 
                     const Rcpp::List tpms,
                     const arma::vec num_cells, 
                     const arma::vec inside, 
                     const double dx, 
                     const arma::vec dt, 
                     const arma::mat sd,
                     const int num_states,
                     const int minstate, 
                     const int maxstate, 
                     const arma::vec entry,
                     arma::field<arma::cube>& lbeta) : n(n), J(J), pr0(pr0), pr_capture(pr_capture), tpms(tpms), num_cells(num_cells), inside(inside), dx(dx), dt(dt), sd(sd), num_states(num_states), minstate(minstate), maxstate(maxstate), entry(entry), lbeta(lbeta) {
    if (num_states > 1) {
      tpm.resize(J); 
      for (int j = 0; j < J - 1; ++j) tpm[j] = Rcpp::as<arma::mat>(tpms[j]); 
    }
    trm.resize(J * num_states); 
    for (int j = 0; j < J - 1; ++j) {
      for (int g = minstate; g < minstate + num_states; ++g) {
        if (sd(j, g - minstate) < 0) continue; 
        trm[g - minstate + j * num_states] = CalcTrm(num_cells, sd(j, g - minstate), dx, inside); 
        trm[g - minstate + j * num_states].t(); 
      }
    }
    pr_cap.resize(n);
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector pr_capvec(pr_capture[i]);
      arma::cube pr_icap(pr_capvec.begin(), num_cells(0), num_states, J, false);
      pr_cap[i] = pr_icap;
    }
  }
  
  void operator()(std::size_t begin, std::size_t end) { 
    
    for (int i = begin; i < end; ++i) {
      arma::mat pr(num_cells(0), num_states, arma::fill::ones);
      pr /= (1.0 * num_cells(0) * num_states); 
      double llk = log((1.0 * num_cells(0) * num_states)); 
      double sum_pr; 
      arma::cube prcap; 
      lbeta(i).slice(J - 1).zeros(); 
      for (int j =  J - 2; j > entry(i) - 1; --j) {
        pr %= pr_cap[i].slice(j + 1);
        if (num_states > 1) {
          pr = tpm[j] * pr.t(); 
        }
        for (int g = minstate; g < minstate + num_states; ++g) {
          if (sd(j, g - minstate) < 0) continue; 
          try {
            pr.col(g) = ExpG(pr.col(g), trm[g - minstate + j * num_states], dt(j));
          } catch(...) {
            break;
          }
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
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_cells number of cells in x,y,total 
//' @param inside 0 if meshpt outside survey region, 1 otherwise 
//' @param dx mesh spacing 
//' @param dt time between occasions 
//' @param sd movement parameter for each occasion 
//' @param num_states 2 = CJS model, 3 = JS model 
//' @param entry time each individual entered survey 
//'
//' @return backwards probabiltiies individual x occasion x mesh x state 
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_movebeta(const int n, const int J,
                                        const arma::mat pr0, 
                                        const Rcpp::List pr_capture, 
                                        const Rcpp::List tpms,
                                        const arma::vec num_cells, 
                                        const arma::vec inside, 
                                        const double dx, 
                                        const arma::vec dt, 
                                        const arma::mat sd, 
                                        const int num_states,
                                        const int minstate, 
                                        const int maxstate, 
                                        const arma::vec entry) {
  
  arma::field<arma::cube> lbeta(n);
  for (int i = 0; i < n; ++i) lbeta(i) = arma::zeros<arma::cube>(num_cells(0), num_states, J); 
  BetaMoveCalculator beta_calculator(n, J, pr0, pr_capture, tpms, num_cells, inside, dx, dt, sd, num_states, minstate, maxstate, entry, lbeta); 
  parallelFor(0, n, beta_calculator, 1); 
  return(lbeta); 
}
