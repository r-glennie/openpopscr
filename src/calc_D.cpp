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
// moveds.cpp: C++ functions to compute abundance at each occasion
//
////////////////////////////////////////////////////////////////////////////////
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>

//' Computes mean density for each occasion
//'
//' @param D density parameter 
//' @param J total number of occasions
//' @param pr0 initial distribution over life states
//' @param tpms list of transition probability matrices (one per occasion)
//'
//' @return  vector with j^th entry the mean density on occasion j 
//'
// [[Rcpp::export]]
arma::vec C_calc_D(  const double D,
                     const int J, 
                   arma::rowvec pr0, 
                   Rcpp::List tpms) {
  
  arma::vec Dt(J); 
  Dt(0) = D * pr0(1); 
  arma::rowvec pr(pr0);
  arma::mat tpm; 
  for (int j = 0; j < J - 1; ++j) {
    tpm = Rcpp::as<arma::mat>(tpms[j]); 
    pr *= tpm; 
    Dt(j + 1) = D * pr(1); 
  }
  return(Dt); 
}
