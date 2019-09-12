# Copyright (c) 2018 Richard Glennie
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################

#' Cormack-Jolly-Seber model class 
#' 
#' @description Cormack-Jolly-Seber model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item scr_data: a ScrData object 
#'   \item start: a named list of starting values 
#'   \item print: (default TRUE) if TRUE then useful output is printed
#' }
#' 
#' Methods include those supplied by ScrModel with the following overwritten: 
#' \itemize{
#'  \item get_par(name, j, k): returns value of parameter "name" for detector j 
#'   on occasion k (if j, k omitted, then returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses as returned by par()
#'  \item calc_initial_distribution(): computes initial distribution over life states (alive, dead)
#'  \item calc_tpms(): returns list of transition probability matrix for each occasion 
#'  \item calc_initial_pdet(): returns probability of being detected somewhere on first occasion seen
#'        for each individual 
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item entry(): return occasion each individual is first detected 
#'  \item fit: fit the model by obtaining the maximum likelihood estimates
#'  \item estimates(): return estimates in a easy to extract list 
#' }
#' 
CjsModel <- R6Class("CjsModel", 
                    inherit = ScrModel, 
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, statemod = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      private$data_ <- data
			index <- 1:data$n_occasions("all")
			if (print) cat("Computing entry occasions for each individual.......")
		  private$entry_ <- apply(data$capthist(), 1, function(x) {min(index[rowSums(x) > 0])}) - 1
		  if (private$data_$n_primary() > 1) {
		    private$entry_ <- private$data_$primary()[private$entry_ + 1] - 1 
		  }
			if (print) cat("done\n")
			if (print) cat("Reading formulae.......")
		  order <- c("phi")
		  private$read_formula(form, detectfn, statemod, order)
		  # add parameters other than detection 
		  private$par_type_[private$detfn_$npars() + 1] <- "p1ms"
		  names(private$form_) <- c(private$detfn_$pars(), "phi")
      # make parameter list 
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi")
      if (print) cat("done\n") 
      if (print) cat("Initilising parameters.......")
      private$initialise_par(start)
      private$read_states() 
      if (print) cat("done\n")
      private$print_ = print
    },
    
    calc_D_llk = function() {warning("No D parameter in CJS model.")}, 
    calc_pdet = function() {warning("No pdet parameter in CJS model.")}, 
  
    calc_initial_pdet = function(pr_capture) {
      inipdet <- rep(0, private$data_$n())
      for (i in 1:private$data_$n()) {
        inipdet[i] <- mean(pr_capture[[i]][,1,private$entry_[i] + 1])
      }
      return(inipdet)
    }, 
    
    calc_initial_distribution = function() {
      n_mesh <- private$data_$n_meshpts()
      nstates <- private$state_$nstates()
      delta <- private$state_$delta() 
      pr0 <- matrix(c(delta, 0), nrow = n_mesh, ncol = nstates + 1, byrow = TRUE)
      pr0 <- pr0 / n_mesh
      return(pr0)
    },
    
    calc_tpms = function() {
      # compute entry probabilities 
      n_occasions <- private$data_$n_occasions()
      n_primary <- private$data_$n_primary() 
      nstates <- self$state()$nstates()
      tpms <- vector("list", length = n_occasions - 1)
      dt <- diff(private$data_$time())
      for (k in 1:(n_occasions - 1)) {
        Q <- matrix(0, nr = nstates + 1, nc = nstates + 1)
        Q[-(nstates+1), -(nstates+1)] <- self$state()$trm(k)
        for (s in 1:nstates) {
          psi <- -log(self$get_par("phi", k = k, m = 1, s = s))
          diag(Q)[s] <- diag(Q)[s] - psi
          Q[s, nstates + 1] <- psi
        }
        tpms[[k]] <- expm(Q * dt[k])
      }
      return(tpms)
    }, 
    
    calc_pr_capture = function() {
      n_primary <- private$data_$n_primary()
      n_occasions <- private$data_$n_occasions("all")
      nstates <- self$state()$nstates()
      kstates <- private$known_states_
      S <- private$data_$n_secondary() 
      if (n_primary == 1) {
        n_primary <- n_occasions
        S <- rep(1, n_occasions)
      }
      enc_rate0 <- self$calc_encrate(transpose = TRUE)
      trap_usage <- usage(private$data_$traps())
      n <- private$data_$n()
      n_meshpts <- private$data_$n_meshpts() 
      n_traps <- private$data_$n_traps()
      capthist <- private$data_$capthist()
      prob <- C_calc_pr_capture(n, 
                                n_occasions, 
                                n_traps, 
                                n_meshpts, 
                                capthist, 
                                enc_rate0, 
                                trap_usage, 
                                nstates,
                                0, 
                                1, 
                                kstates, 
                                self$data()$detector_type(), 
                                n_primary, 
                                S,
                                private$entry_)
      return(prob)
    },
    
    calc_llk = function(param = NULL, names = NULL) {
      if (!is.null(names)) names(param) <- names 
      if (!is.null(param)) {
        slen <- length(self$state()$par())
        param2 <- param 
        if (slen > 0) {
          ind <- seq(length(param) - slen + 1, length(param))
          self$state()$set_par(param[ind])
          param2 <- param[-ind]
        }
        self$set_par(private$convert_vec2par(param2));
      }
      # compute transition probability matrices 
      tpms <- self$calc_tpms()
      # initial distribution 
      pr0 <- self$calc_initial_distribution()
      # compute probability of capture histories 
      # across all individuals, occasions and traps 
      pr_capture <- self$calc_pr_capture()
      # compute likelihood for each individual
      n <- private$data_$n()
      n_occasions <- private$data_$n_occasions()
      n_meshpts <- private$data_$n_meshpts() 
      nstates <- self$state()$nstates() + 1
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, private$entry_)
      # compute probability of initial detection
      inipdet <- self$calc_initial_pdet(pr_capture) 
      llk <- llk - sum(log(inipdet))
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
  
  entry = function() {return(private$entry_)}, 
  
  estimates = function() {
      ests <- NULL
      if (is.null(private$mle_)) {
        ests <- "Fit model using $fit method"
      } else {
        ests$par <- private$results_
      }
      return(ests)
    }
),
                   
  private = list(
		entry_ = NULL,
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                        list(start[[names[i]]]))
      }
      private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
      # compute initial parameters for each jkm
      private$compute_par()
      return(invisible())
    }, 
		
		calc_forwback = function(forw = NULL, back = NULL) {
		  if (is.null(forw) & is.null(back)) forw <- back <- TRUE
		  if (is.null(forw)) forw <- FALSE
		  if (is.null(back)) back <- FALSE
		  # initial distribution 
		  pr0 <- self$calc_initial_distribution()
		  # compute probability of capture histories 
		  # across all individuals, occasions and traps 
		  pr_capture <- self$calc_pr_capture()
		  # compute lalpha for each individual
		  n <- private$data_$n()
		  n_occasions <- private$data_$n_occasions()
		  n_meshpts <- private$data_$n_meshpts() 
		  # get tpms for state model 
		  nstates <- self$state()$nstates() 
		  tpms <- self$calc_tpms()
		  # compute forward-backward 
		  if (forw) lalpha <- C_calc_alpha(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, private$entry_)
		  if (back) lbeta <- C_calc_beta(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, private$entry_) 
		  if (forw & back) return(list(lalpha = lalpha, lbeta = lbeta))
		  if (forw) return(lalpha)
		  if (back) return(lbeta)
		  return(0)
		}, 
		
		read_states = function() {
		  nstates <- self$state()$nstates() + 1 
		  kstates <- array(1, dim = c(private$data_$n(), private$data_$n_occasions("all"), nstates))
		  covtypes <- private$data_$get_cov_list()$cov_type
		  snms <- self$state()$names()
		  grpnms <- self$state()$groupnms()
		  if ("dead" %in% grpnms) stop("Cannot have a state variable named 'dead'. This is a reserved word.")
		  grps <- self$state()$groups()
		  if ("i" %in% covtypes |  "ik" %in% covtypes) {
		    wh <- min(which(covtypes %in% c("i", "ik")))
		    cov <- private$data_$get_cov_list()$cov[[wh]]
		    type <- covtypes[wh]
		    for (i in 1:private$data_$n()) {
		      for (k in 1:private$data_$n_occasions()) {
		        s <- private$data_$covs(i = i, k = k)
		        for (g in 1:length(grpnms)) {
		          if (grpnms[g] %in% names(s)) {
		            occu <- grps[,g] %in% s[[grpnms[[g]]]]
		            if (any(occu)) kstates[i, k, !occu] <- -1
		          }
		        }
		      }
		      if ("dead" %in% names(s)) {
		        if(!is.na(s$dead)) kstates[i, k, -(nstates + 1)] <- -1
		      }
		    }
		  }
		  private$known_states_ <- kstates
		  invisible()
		}
  )                 
)



