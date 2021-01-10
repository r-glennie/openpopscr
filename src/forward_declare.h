// Copyright (c) 2018 Richard Glennie, University of St Andrews
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
#ifndef OPENPOPSCR_FORWARD_DECALRE_H
#define OPENPOPSCR_FORWARD_DECALRE_H

#include <iostream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>

using namespace RcppParallel; 

arma::sp_mat CalcTrm(const arma::vec num_cells, const double sd, const double dx, const arma::mat inside); 
arma::vec ExpG(const arma::vec& v_in,
               const arma::sp_mat& a,
               const double& t,
               const int& krylov_dim = 30,
               const double& tol = 1e-10); 

#endif 