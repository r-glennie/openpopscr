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
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <iostream>
#include <RcppArmadillo.h>

arma::sp_mat CalcTrm(const arma::vec num_cells, const double sd, const double dx, const arma::vec inside) {
  arma::sp_mat tpr = arma::zeros<arma::sp_mat>(num_cells(0), num_cells(0));
  double rate = sd * sd / (2 * dx * dx);
  int s;
  double sum; 
  for (int i = 0; i < num_cells(1); ++i) {
    for (int j = 0; j < num_cells(2); ++j) {
      s = i + num_cells(1) * j; 
      sum = 0; 
      if (inside(s) > 0) {
        if (i < num_cells(1) - 1) {
          if (inside(s + 1) > 0) {
            tpr(s, s + 1) = rate; 
            sum += rate; 
          }
        }
        if (i > 0) {
          if (inside(s - 1) > 0) {
            tpr(s, s - 1) = rate; 
            sum += rate; 
          }
        }
        if (j < num_cells(2) - 1) {
          if (inside(s + num_cells(1)) > 0) {
            tpr(s, s + num_cells(1)) = rate; 
            sum += rate; 
          }
        }
        if (j > 0) {
          if (inside(s - num_cells(1)) > 0) {
            tpr(s, s - num_cells(1)) = rate; 
            sum += rate;
          }
        }
      }
      tpr(s, s) = -sum; 
    }
  }
  return tpr.t();
}

arma::vec ExpG(const arma::vec v,
               const arma::sp_mat a,
               const double& t,
               const int& krylov_dim = 30,
               const double& tol = 1e-16) {
  double m = fmin(a.n_rows, krylov_dim);
  double anorm = arma::norm(a, "Inf");
  double mxrej = 10;
  double mx;
  double btol = 1e-7;
  double gamma = 0.9;
  double mb = m;
  int nstep = 0;
  double t_now = 0;
  double t_step;
  double delta = 1.2;
  double t_out = fabs(t);
  double s_error = 0;
  double rndoff = anorm * 1e-16;
  
  int k1 = 1;
  double xm = 1 / m;
  double normv = norm(v);
  double avnorm;
  double beta = normv;
  double fact = std::pow((m + 1) / std::exp(1), m + 1) * std::sqrt(2 * M_PI * (m + 1));
  double t_new = (1.0 / anorm) * std::pow((fact * tol) / (4 * beta * anorm), xm);
  double s = std::pow(10, std::floor(std::log10(t_new)) - 1);
  t_new = std::ceil(t_new / s) * s;
  double sgn = t > 0 ? 1 : -1;
  int ireject;
  double err_loc;
  double phi1;
  double phi2;
  
  arma::vec w = v;
  double hump = normv;
  arma::mat vmat = arma::zeros<arma::mat>(a.n_rows, m + 1);
  arma::mat hmat = arma::zeros<arma::mat>(m + 2, m + 2);
  arma::mat fmat;
  arma::vec p;
  while (t_now < t_out) {
    Rcpp::checkUserInterrupt();
    ++nstep;
    t_step = fmin(t_out - t_now, t_new);
    vmat.zeros();
    hmat.zeros();
    vmat.col(0) = (1 / beta) * w;
    for (int j = 0; j < m; ++j) {
      p = a * vmat.col(j);
      for (int i = 0; i <= j; ++i) {
        hmat(i, j) = arma::dot(vmat.col(i), p);
        p -= hmat(i, j) * vmat.col(i);
      }
      s = norm(p);
      if (s < btol) {
        k1 = 0;
        mb = j;
        t_step = t_out - t_now;
        break;
      }
      hmat(j + 1, j) = s;
      vmat.col(j + 1) = (1 / s) * p;
    }
    if (k1 != 0) {
      hmat(m + 1, m) = 1;
      avnorm = arma::norm(a * vmat.col(m));
    }
    ireject = 0;
    while (ireject <= mxrej) {
      mx = mb + k1;
      fmat = arma::expmat(sgn * t_step * hmat.submat(0, 0, mx, mx));
      if (k1 == 0) {
        err_loc = btol;
        break;
      }
      else {
        phi1 = fabs(beta * fmat(m, 0));
        phi2 = fabs(beta * fmat(m + 1, 0) * avnorm);
        if (phi1 > 10 * phi2) {
          err_loc = phi2;
          xm = 1 / m;
        }
        else if (phi1 > phi2) {
          err_loc = (phi1 * phi2) / (phi1 - phi2);
          xm = 1 / m;
        }
        else {
          err_loc = phi1;
          xm = 1 / (m - 1);
        }
      }
      if (err_loc <= delta * t_step * tol) break;
      else {
        t_step = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
        s = std::pow(10, std::floor(std::log10(t_step)) - 1);
        t_step = std::ceil(t_step / s) * s;
        if (ireject == mxrej) {
          //std::cout << "error: requested tolerance too high for Krylov approximation" << std::endl;
        }
        ++ireject;
      }
    }
    mx = mb + fmax(0, k1 - 1);
    w = vmat.cols(0, mx) * beta * fmat.col(0).rows(0, mx);
    beta = arma::norm(w);
    hump = fmax(hump, beta);
    
    t_now = t_now + t_step;
    t_new = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
    s = std::pow(10, std::floor(std::log10(t_new) - 1));
    t_new = std::ceil(t_new / s) * s;
    
    err_loc = fmax(err_loc, rndoff);
    s_error += err_loc;
  }
  double err = s_error;
  hump = hump / normv;
  return w;
}

//' Computes log-likelihood of Jolly-Seber model 
//'
//' @param n number of individuals 
//' @param J total number of occasions 
//' @param M total number of mesh points
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_cores number of processor cores to use in parallelisation 
//' @param num_states 2 = CJS model, 3 = JS model 
//'
//' @return log-likelihood value 
//' 
// [[Rcpp::export]]
double C_calc_move_llk(const int n, const int J,
                  const arma::mat pr0, 
                  const Rcpp::List pr_capture, 
                  const Rcpp::List tpms,
                  const arma::vec num_cells, 
                  const arma::vec inside, 
                  const double dx, 
                  const arma::vec dt, 
                  const arma::vec sd, 
                  const int num_cores,
                  const int num_states,
                  const arma::vec entry) {
  
  arma::vec illk(n);
  arma::mat pr; 
  double llk; 
  double sum_pr; 
  arma::mat tpm; 
  Rcpp::NumericVector pr_capvec; 
  arma::sp_mat trm(num_cells(0), num_cells(0));  
  int alive_col = 0; 
  if (num_states > 2) alive_col = 1; 
  int M = num_cells(0); 
  for (int i = 0; i < n; ++i) {
    llk = 0; 
    sum_pr; 
    pr = pr0; 
    tpm; 
    pr_capvec = pr_capture[i]; 
    arma::cube pr_cap(pr_capvec.begin(), M, num_states, J, false); 
    for (int j = entry(i); j < J - 1; ++j) {
      pr %= pr_cap.slice(j); 
      trm = CalcTrm(num_cells, sd(j), dx, inside);
      if (num_states > 1) {
        tpm = Rcpp::as<arma::mat>(tpms[j]); 
        pr *= tpm; 
      }
      try {
        pr.col(alive_col) = ExpG(pr.col(alive_col), trm, dt(j)); 
      } catch(...) {
        return -arma::datum::inf; 
      }
      sum_pr = accu(pr); 
      llk += log(sum_pr); 
      pr /= sum_pr; 
    }
    pr %= pr_cap.slice(J - 1);
    llk += log(accu(pr)); 
    illk(i) = llk;  
  }
  return(arma::accu(illk)); 
}

//' Computes detection probability (seen at least once) for Jolly-Seber model 
//'
//' @param J total number of occasions 
//' @param pr0 initial distribution over life states
//' @param pr_captures list of empty capture histories, see calc_pdet() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//'
//' @return pdet = probability seen at some time on the survey 
//' 
// [[Rcpp::export]]
double C_calc_move_pdet(const int J, 
                   arma::mat pr0, 
                   Rcpp::List pr_captures,
                   Rcpp::List tpms,
                   const arma::vec num_cells, 
                   const arma::vec inside, 
                   const double dx, 
                   const arma::vec dt,
                   const arma::vec sd, 
                   const int num_states) {
  
  double pdet = 0; 
  arma::mat pr(pr0);
  arma::sp_mat trm(num_cells(0), num_cells(0)); 
  double sum_pr;
  arma::mat tpm;
  arma::mat pr_capture;
  int alive_col = 0; 
  if (num_states > 2) alive_col = 1; 
  for (int j = 0; j < J - 1; ++j) {
    pr_capture = Rcpp::as<arma::mat>(pr_captures[j]);
    pr %= pr_capture;
    trm = CalcTrm(num_cells, sd(j), dx, inside);
    if (num_states > 1) {
      tpm = Rcpp::as<arma::mat>(tpms[j]); 
      pr *= tpm; 
    }
    try {
      pr.col(alive_col) = ExpG(pr.col(alive_col), trm, dt(j)); 
    } catch(...) {
      return -arma::datum::inf; 
    }
    sum_pr = accu(pr); 
    pdet += log(sum_pr); 
    pr /= sum_pr; 
  }
  pr %= pr_capture;
  pdet += log(accu(pr)); 
  pdet = 1 - exp(pdet); 
  return(pdet); 
}



