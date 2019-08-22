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
  const int J; // number of occasions 
  const int K; // number of traps 
  const int M; // number of mesh points 
  const int alive_col; // column that contains alive state 
  const arma::cube& capthist; // capthist history records: individuals x occasion x trap 
  const Rcpp::List enc; // encounter rate: occasion x mesh point x trap
  const arma::mat& usage; // usage of traps: trap x occasion 
  const int num_states; // number of hidden states in life history model 
  const arma::cube& known_state; // -1 if known not to be in that state 
  const int minstate; // states before alive state 
  const int maxstate; // states after alives states
  const int detector_type; // type of detector, see initialize in ScrData  
  const int n_prim; // number of primary occasions 
  const arma::vec& S; // number of secondary occasion per primary 
  const arma::vec& entry; // occasion each individual enters survey 
  
  // derived variables 
  std::vector<arma::cube> enc0; 
  std::vector<arma::mat> total_enc; 
  std::vector<arma::mat> log_total_enc; 
  std::vector<arma::mat> log_total_penc; 
  std::vector<arma::cube> logenc0; 
  std::vector<arma::cube> log_penc; 

  // output 
  // capture probability for record of individual x occasion x mesh point
  arma::field<arma::cube>& probfield;
  
  // initialise
  PrCaptureCalculator(const int J, 
                        const int K, 
                        const int M, 
                        const int alive_col, 
                        const arma::cube& capthist, 
                        const Rcpp::List enc, 
                        const arma::mat& usage, 
                        const int num_states,
                        const int minstate, 
                        const int maxstate, 
                        const arma::cube& known_state, 
                        const int detector_type,
                        const int n_prim, 
                        const arma::vec& S, 
                        const arma::vec& entry, 
                        arma::field<arma::cube>& probfield) : J(J), K(K), M(M), 
                        alive_col(alive_col), 
                        capthist(capthist), 
                        enc(enc), 
                        usage(usage), 
                        num_states(num_states), 
                        minstate(minstate), 
                        maxstate(maxstate), 
                        known_state(known_state), 
                        detector_type(detector_type),
                        n_prim(n_prim), 
                        S(S), 
                        entry(entry), 
                        probfield(probfield) {
    
    //// compute dervied quantities 
    // encounter rate per state 
    enc0.resize(num_states); 
    for (int g = 0; g < num_states; ++g) {
      enc0[g] = Rcpp::as<arma::cube>(enc[g]); 
    }
    if (detector_type == 3) {
      total_enc.resize(num_states);
      log_total_enc.resize(num_states);
      log_total_penc.resize(num_states);
    } else if (detector_type == 2) {
      log_penc.resize(num_states);
    }
    for (int g = 0; g < num_states; ++g) {
      if (detector_type == 3) { 
        // multi detector: 
        //  log_total_enc: log total hazard over detectors for each mesh point x occasion 
        //  log_total_penc: log total probability of detection for each mesh point x occasion 
        total_enc[g] = arma::zeros<arma::mat>(M, J); 
        log_total_enc[g] = arma::zeros<arma::mat>(M, J); 
        log_total_penc[g] = arma::zeros<arma::mat>(M, J); 
        int j = -1; 
        for (int prim = 0; prim < n_prim; ++prim) {
          for (int s = 0; s < S(prim); ++s) {
            ++j; 
            total_enc[g].col(j) = enc0[g].slice(j) * usage.col(j); 
          }
        }
        log_total_enc[g] = log(total_enc[g]); 
        log_total_penc[g] = 1.0 - exp(-total_enc[g]); 
        log_total_penc[g] = log(log_total_penc[g] + 1e-16); 
      } else if (detector_type == 2) {
        // proximity detector: 
        //  log_penc: log-probability seen at some point mesh point x trap x occasion 
        log_penc[g] = arma::zeros<arma::cube>(M, K, J); 
        int j = -1; 
        for (int prim = 0; prim < n_prim; ++prim) {
          for (int s = 0; s < S(prim); ++s) {
            ++j; 
            for (int k = 0; k < K; ++k) { 
              log_penc[g].slice(j).col(k) = 1.0 - exp(-enc0[g].slice(j).col(k) * usage(k, j)); 
              log_penc[g].slice(j).col(k) = log(log_penc[g].slice(j).col(k) + 1e-16); 
            }
          }
        }
      }
      // log encounter rate occasion x mesh point x trap
      logenc0.resize(num_states);
      logenc0[g] = log(enc0[g]); 
    }
  } 
  
  void operator()(std::size_t begin, std::size_t end) { 
    // loop over individuals 
    for (int i = begin; i < end; ++i) {
      arma::vec savedenc(M); 
      double sumcap; 
      int g = 0; // current state processed 
      int j = -1; // current occasion processed 
      // loop over primary occasions
      for (int prim = 0; prim < n_prim; ++prim) {
        bool unseen = true; // individual seen or not in this primary?
        if (prim > entry(i) - 1) {
          // loop over secondary occasions 
          for (int s = 0; s < S(prim); ++s) {
            ++j; // increment occasion number 
            // loop over detectable states
            for (int gp = minstate; gp < minstate + num_states; ++gp) {
              g = gp - minstate; // current state 
              // if state known to be impossible, leave probfield zero 
              if (known_state(i, j, gp) < 0) {
                continue; 
              } else {
                // state is possible 
                if (detector_type == 3) savedenc.zeros(); 
                sumcap = 0; 
                // independent detectors 
                if (detector_type != 3) {
                  // loop over detectors 
                  for (int k = 0; k < K; ++k) {
                    if (usage(k, j) < 1e-16) continue; // detector unused 
                    if (detector_type == 1 | detector_type == 4) {
                      // count detector Poisson
                      probfield(i).slice(prim).col(gp) += capthist(i, j, k) * logenc0[g].slice(j).col(k) - usage(k, j) * enc0[g].slice(j).col(k);
                    } else if (detector_type == 2) {
                      // proximity detector 
                      probfield(i).slice(prim).col(gp) += capthist(i, j, k) * log_penc[g].slice(j).col(k) - (1.0 - capthist(i, j, k)) * usage(k, j) * enc0[g].slice(j).col(k);
                    }
                    // seen? 
                    if (capthist(i, j, k) > 1e-16) unseen = false; 
                  }
                }
                // dependent detectors 
                if (detector_type == 3) {
                  // get detector record for this individual and occasion 
                  arma::vec cap_ij = capthist(arma::span(i), arma::span(j), arma::span::all);  
                  sumcap = arma::accu(cap_ij);
                  // seen? 
                  if (sumcap > 0) unseen = false; 
                  savedenc += logenc0[g].slice(j) * cap_ij;  
                  if (!unseen) probfield(i).slice(prim).col(gp) += savedenc - sumcap * log_total_enc[g].col(j); 
                  probfield(i).slice(prim).col(gp) += -(1.0 - sumcap) * total_enc[g].col(j) + sumcap * log_total_penc[g].col(j); 
                }
              }
              probfield(i).slice(prim).col(gp) = exp(probfield(i).slice(prim).col(gp));
            }
          }
        } else {
          // if not entered in this primary yet 
          j = j + S(prim); // move occasion number to end of primary 
        }
        if (unseen) {
          for (int gpi = 0; gpi < minstate; ++gpi) probfield(i).slice(prim).col(gpi).ones(); 
          for (int gpi = minstate + num_states; gpi < minstate + num_states + maxstate; ++gpi) probfield(i).slice(prim).col(gpi).ones(); 
        }
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
//' @param capthist capthist array 
//' @param enc0 encounter rate array, see calc_pr_capture() in JsModel
//' @param usage matrix with J x K where (j,k) entry is usage of trap k in occasion j
//' @param num_states number of alive states 
//' @param minstate number of states before alive (Scr,Cjs = 0, JS = 1)
//' @param maxstate number of states after alive (Scr = 0, Cjs/js = 1)
//' @param detector_type 1 = count, 2 = proximity/binary, 3 = multi-catch, 4 = transect 
//' @param n_prim number of primary occasions 
//' @param S number of secondary occasions per primary occasion 
//' @param entry occasion each individual entered survey 
//'
//' @return  Array with (i,j,m) entry the probability of capture record for individual i in occasion j given activity centre at mesh point m  
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_pr_capture(const int n, const int J, const int K, const int M, 
                             const arma::cube& capthist, 
                             const Rcpp::List enc0,
                             const arma::mat usage, 
                             const int num_states,
                             const int minstate, 
                             const int maxstate,
                             const arma::cube& known_state, 
                             const int detector_type, 
                             const int n_prim, 
                             const arma::vec S, 
                             const arma::vec entry) {
  int alive_col = 1; 
  if (num_states < 3) alive_col = 0; 
  arma::field<arma::cube> probfield(n);
  for (int i = 0; i < n; ++i) probfield(i) = arma::zeros<arma::cube>(M, num_states, n_prim); 
  PrCaptureCalculator pr_capture_calc(J, K, M, alive_col, capthist, enc0, usage, num_states, minstate, maxstate, known_state, detector_type, n_prim, S, entry, probfield); 
  parallelFor(0, n, pr_capture_calc, 1); 
  return(probfield);
}
