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

#' SCR model class 
#' 
#' @description Spatial capture-recapture model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item scr_data: a ScrData object 
#'   \item start: a named list of starting values 
#'   \item num_cores (optional, default = 1): number of processors cores to use 
#'   in parallelised code 
#'   \item print (default = TRUE): if TRUE then helpful output is printed
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k (if j, k, then returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses. Note, the model will simulate 
#'    data using this parameter, but will only present inference based on the maximum likelihood
#'    estimates. 
#'  \item set_mle(mle, V, llk): set maximum likelihood point with parameters mle, 
#'  covariance matrix V, and value llk
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit: fit the model by obtaining the maximum likelihood estimates. Estimates of
#'        density are obtained from parametric boostrap with nsim resamples. 
#'  \item simulate(): simulate ScrData object from fitted model
#'  \item par(): return current parameter of the model 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#' }
#' 
ScrModel <- R6Class("ScrModel", 
  public = list(
    
    initialize = function(form, data, start, num_cores = 1, print = TRUE) {
      private$data_ <- data
      if (print) cat("Reading formulae.......")
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$form_[[1]]<- form[par_names == "lambda0"][[1]]
      private$form_[[2]] <- form[par_names == "sigma"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c("lambda0", "sigma")
      private$make_par() 
      private$link2response_ <- list(lambda0 = "exp", 
                            sigma = "exp", 
                            D = "exp")
      private$response2link_ <- list(lambda0 = "log", 
                                     sigma = "log", 
                                     D = "log")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$num_cores_ = num_cores
      private$print_ = print 
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
      if (name == "D") return(do.call(private$link2response_$D, list(as.numeric(self$par()["D"]))))
     covs <- private$data_$covs(j = j, k = k, m = m)
     i_par <- which(names(private$form_) == name) 
     if (is.null(i_par)) stop("No parameters with that name.")
     X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
     theta <- private$par_[[i_par]]
     l_par <- which(names(private$link2response_) == name)
     resp <- do.call(private$link2response_[[l_par]], list(X %*% theta))
     return(resp)
    }, 
    
    set_par = function(par) {
      private$par_ <- par
    },
    
    set_mle = function(mle, V, llk) {
      private$mle_ <- mle
      private$par_ <- mle
      private$llk_ <- -llk
      private$V_ <- V
      private$make_results()
    }, 
    
    calc_D_llk = function() {
      D <- do.call(private$link2response_$D, list(self$par()$D))
      A <- private$data_$area()
      n <- private$data_$n()
      pdet <- self$calc_pdet()
      llk <- n * log(D * A * pdet) - D * A * pdet - lfactorial(n)
      names(llk) <- NULL 
      return(llk)
    },
    
    calc_initial_distribution = function() {
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(1, nrow = n_mesh, ncol = 1, byrow = TRUE)
      pr0 <- pr0 / n_mesh
      return(pr0)
    },
    
    calc_pr_capture = function() {
      dist <- t(private$data_$distances())
      n_occasions <- private$data_$n_occasions()
      enc_rate0 <- array(0, dim = c(nrow(dist), ncol(dist), n_occasions)) 
      for (k in 1:n_occasions) {
        lambda0 <- as.vector(self$get_par("lambda0", k = k, m = 1))
        sigma <- as.vector(self$get_par("sigma", k = k, m = 1))
        enc_rate0[,,k] <- lambda0 * exp(-dist ^ 2 / (2  * sigma ^ 2))
      }
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
                                private$num_cores_, 
                                1, 
                                self$data()$detector_type(), 
                                n_occasions, 
                                rep(1, n_occasions),
                                rep(0, n))
      return(prob)
    },
    
    calc_pdet = function() {
      # compute probability of zero capture history 
      dist <- private$data_$distances()
      n_occasions <- private$data_$n_occasions()
      enc_rate <- array(0, dim = c(n_occasions, nrow(dist), ncol(dist))) 
      for (k in 1:n_occasions) {
        lambda0 <- as.vector(self$get_par("lambda0", k = k, m = 1))
        sigma <- as.vector(self$get_par("sigma", k = k, m = 1))
        enc_rate[k,,] <- lambda0 * exp(-dist ^ 2 / (2  * sigma ^ 2))
      }
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      for (j in 1:private$data_$n_occasions()) {
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = 1)
        pr_empty[[j]][, 1] <- exp(-t(trap_usage[, j]) %*% enc_rate[j,,])
      }
      pr0 <- self$calc_initial_distribution()
      tpms <- list(matrix(0, nr = 2, nc = 2))
      pdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms, 1); 
      return(pdet)
    },
    
    calc_llk = function(param = NULL, names = NULL) {
      if (!is.null(names)) names(param) <- names 
      if (!is.null(param)) self$set_par(private$convert_vec2par(param));
      # initial distribution 
      pr0 <- self$calc_initial_distribution()
      # compute probability of capture histories 
      # across all individuals, occasions and traps 
      pr_capture <- self$calc_pr_capture()
      # compute likelihood for each individual
      n <- private$data_$n()
      n_occasions <- private$data_$n_occasions()
      n_meshpts <- private$data_$n_meshpts() 
      tpms <- list(matrix(0, nr = 2, nc = 2))
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms,
			private$num_cores_, 1, rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_pdet())
      llk <- llk + self$calc_D_llk()
      cat("llk:", llk, "\n")
      return(llk)
    },
    
    fit = function(ini_par = NULL, nlm.args = NULL) {
      if (!is.null(ini_par)) self$set_par(ini_par)
      par <- self$par()
      w_par <- private$convert_par2vec(par)
      t0 <- Sys.time()
      if (private$print_) cat("Fitting model..........\n")
      args <- c(list(private$calc_negllk, w_par, names = names(w_par), hessian = TRUE), nlm.args)
      mod <- do.call(nlm, args)
      t1 <- Sys.time()
      difft <- t1 - t0 
      if (private$print_) cat("Completed model fitting in", difft, attr(difft, "units"), "\n")
      code <- mod$code
      if (code > 2) warning("model failed to converge with nlm code ", code)
      if (private$print_ & code < 3) cat("Checking convergence.......converged", "\n")
      mle <- mod$estimate
      names(mle) <- names(w_par)
      mle <- private$convert_vec2par(mle)
      #mle <- lapply(mle, function(x) {y <- x; names(y) <- NULL; return(y)})
      self$set_par(mle)
      private$mle_ <- mle
      private$llk_ <- -mod$minimum
      if (private$print_) cat("Computing variance.......")
      private$V_ <- solve(mod$hessian)
      if (private$print_) cat("done\n")
      if (any(diag(private$V) <= 0)) {
        cat("Variance estimates not reliable, do a bootstrap.")
        return(0)
      } else {
        sd <- sqrt(diag(private$V_))
        names(sd) <- names(w_par)
        private$make_results()
      }
    }, 
    
  print = function() {
    options(scipen = 999)
    if (is.null(private$mle_)) {
      print("Fit model using $fit method")
    } else {
      cat("PARAMETER ESTIMATES (link scale)\n")
      print(signif(private$results_, 4))
    }
    options(scipen = 0)
  }, 
  
  simulate = function(seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    new_dat <- simulate_scr(self$par(), 
                            self$data()$n_occasions(), 
                            self$data()$traps(), 
                            self$data()$mesh(), 
                            self$data()$time(), 
                            seed, 
                            private$print_)
    return(new_dat)
  }, 
  
  par = function() {return(private$par_)},
  mle = function() {return(private$mle_)},
  data = function() {return(private$data_)}, 
  
  estimates = function() {
      ests <- NULL
      if (is.null(private$mle_)) {
        ests <- "Fit model using $fit method"
      } else {
        ests$par <- private$results_
        ests$D <- private$D_tab_
      }
      return(ests)
    },
  
    cov_matrix = function() {return(private$V_)}, 
    mle_llk = function() {return(private$llk_)}
  
),
                   
  private = list(
    data_ = NULL,
    form_ = NULL, 
    par_ = NULL, 
    link2response_ = NULL, 
    response2link_ = NULL, 
    mle_ = NULL,
    results_ = NULL,
    D_tab_ = NULL, 
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05, 
    num_cores_ = NULL,
    print_  = NULL, 
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_par <- numeric(4)
      private$par_ <- vector(mode = "list", length = 3)
      for (par in 1:2) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        par_vec <- rep(0, n_par[par])
        names(par_vec) <- colnames(X)
        private$par_[[par]] <- par_vec
      }
      private$par_[[3]] <- 0
      names(private$par_) <- c(names(private$form_), "D") 
    }, 
    
    initialise_par = function(start) {
        private$par_$lambda0[1] <- do.call(private$response2link_$lambda0, 
                                           list(start$lambda0))
        private$par_$sigma[1] <-do.call(private$response2link_$sigma, 
                                           list(start$sigma))
        private$par_$D <- do.call(private$response2link_$D, 
                                           list(start$D))
    }, 
    
    make_results = function() {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      par <- private$convert_par2vec(private$mle_)
      sd <- sqrt(diag(private$V_))
      if (private$print_) cat("Computing confidence intervals.......")
      ci <- private$calc_confint()
      if (private$print_) cat("done\n")
      results <- cbind(par, sd, ci$LCL, ci$UCL)
      colnames(results) <- c("Estimate", "Std. Error", "LCL", "UCL")
      rownames(results) <- names(par)
      private$results_ <- results 
    
  }, 
  
   calc_confint = function() {
      V <- private$V_ 
      sds <- sqrt(diag(V))
      est <- private$convert_par2vec(private$mle_)
      alp <- qnorm(1 - private$sig_level_ / 2)
      lcl <- est - alp * sds
      ucl <- est + alp * sds 
      names(lcl) <- names(ucl) <- names(est)
      return(list(LCL = lcl, UCL = ucl))
    },
  
     calc_negllk = function(param = NULL, names = NULL) {
       negllk <- -self$calc_llk(param, names)
       return(negllk)
    },
  
   convert_par2vec = function(par) {
      return(unlist(par))
    },
    
    convert_vec2par = function(vec) {
      par <- NULL
      n_occasions <- private$data_$n_occasions()
      names <- names(vec)
      par$lambda0 <- vec[grep("lambda0", names)]
      names(par$lambda0) <- gsub("lambda0.", "", names(par$lambda0))
      par$sigma <- vec[grep("sigma", names)]
      names(par$sigma) <- gsub("sigma.", "", names(par$sigma))
      par$D <- vec["D"]
      names(par$D) <- NULL 
      return(par)
    }
  )                 
)



