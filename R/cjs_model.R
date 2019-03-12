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
#'   \item print: (defualt TRUE) if TRUE then useful output is printed
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k): returns value of parameter "name" for detector j 
#'   on occasion k (if j, k omitted, then returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses. Note, the model will simulate 
#'    data using this parameter, but will only present inference based on the maximum likelihood
#'    estimates. 
#'  \item set_mle(mle, V, llk): sets model at maximum likelihood with parameters mle, 
#'  covariance matrix V, and likelihood value llk 
#'  \item calc_initial_distribution(): computes initial distribution over life states (alive, dead)
#'  \item calc_tpms(): returns list of transition probability matrix for each occasion 
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit: fit the model by obtaining the maximum likelihood estimates
#'  \item simulate(): simulate ScrData object from fitted model
#'  \item par(): return current parameter of the model 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#' }
#' 
CjsModel <- R6Class("CjsModel", 
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, print = TRUE) {
      private$data_ <- data
			index <- 1:data$n_occasions("all")
			if (print) cat("Computing entry occasions for each individual.......")
		  private$entry_ <- apply(data$capthist(), 1, function(x) {min(index[rowSums(x) > 0])}) - 1
		  if (private$data_$n_primary() > 1) {
		    private$entry_ <- private$data_$primary()[private$entry_ + 1] - 1 
		  }
			if (print) cat("done\n")
			if (print) cat("Reading formulae.......")
		  private$detectfn_ <- detectfn 
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      # detection function 
      if (is.null(detectfn)) {
        private$detfn_ <- DetFn$new()
      } else if (class(detectfn)[1] == "character") {
        private$detfn_ <- DetFn$new(fn = detectfn)
      } else {
        private$detfn_ <- detectfn 
      }
      for (i in 1:private$detfn_$npars()) {
        find <- par_names == private$detfn_$par(i)
        if (all(!find)) stop("Parameters in formulae incorrect.")
        private$form_[[i]]<- form[find][[1]]
      }
      private$form_[[private$detfn_$npars() + 1]] <- form[par_names == "phi"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c(private$detfn_$pars(), "phi")
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi")
      if (print) cat("done\n") 
      if (print) cat("Initilising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$print_ = print
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
     if (name == "phi" & is.null(k)) {
       k <- 1:(private$data_$n_occasions() - 1) 
       n_primary <- private$data_$n_primary() 
       if (n_primary > 1) k <- match(2:n_primary, private$data_$primary())
     }
     covs <- private$data_$covs(j = j, k = k, m = m)
     i_par <- which(names(private$form_) == name) 
     if (length(i_par) == 0) stop("No parameter with that name.")
     X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
     if (name %in% c("phi") & "t" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("t1", sep ="")] 
     }
     if (name %in% c("phi") & "primary" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("primary1", sep ="")] 
     }
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
      private$llk_ <- -llk
      private$V_ <- V
      private$make_results()
    }, 
    
    calc_initial_distribution = function() {
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1, 0), nrow = n_mesh, ncol = 2, byrow = TRUE)
      pr0 <- pr0 / n_mesh
      return(pr0)
    },
    
    calc_tpms = function() {
      # compute entry probabilities 
      n_occasions <- private$data_$n_occasions()
      n_primary <- private$data_$n_primary() 
      if (n_primary > 1) first <- match(2:n_primary, private$data_$primary())
      tpms <- vector("list", length = n_occasions)
      dt <- diff(private$data_$time())
      for (k in 1:(n_occasions - 1)) {
        occ <- k 
        if (n_primary > 1) occ <- first[k]
        phi <- self$get_par("phi", j = 1, k = occ, m = 1)
        phi <- phi^dt[k]
        tpms[[k]] <- matrix(c(phi, 0, 
                            1 - phi, 1), nrow = 2, ncol = 2)
      }
      return(tpms)
    }, 
    
    calc_encrate = function(transpose = FALSE) {
      dist <- private$data_$distances()
      if (transpose) dist <- t(dist)
      n_occasions <- private$data_$n_occasions()
      if(transpose) enc_rate <- array(0, dim = c(nrow(dist), ncol(dist), n_occasions)) 
      if(!transpose) enc_rate <- array(0, dim = c(n_occasions, nrow(dist), ncol(dist))) 
      n_det_par <- self$detectfn()$npars()
      det_par <- vector(mode = "list", length = n_det_par)
      for (k in 1:n_occasions) {
        for (dpar in 1:n_det_par) det_par[[dpar]] <- as.vector(self$get_par(self$detectfn()$par(dpar), k = k, m = 1))
        if(transpose) enc_rate[,,k] <- self$detectfn()$h(dist, det_par)
        if(!transpose) enc_rate[k,,] <- self$detectfn()$h(dist, det_par)
      }
      # add epsilon to stop log(0.0)
      enc_rate <- enc_rate + 1e-16
      return(enc_rate) 
    }, 
    
    calc_pr_capture = function() {
      n_primary <- private$data_$n_primary()
      n_occasions <- private$data_$n_occasions("all")
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
                                2, 
                                self$data()$detector_type(), 
                                n_primary, 
                                S,
                                private$entry_)
      return(prob)
    },
    
    calc_llk = function(param = NULL, names = NULL) {
      if (!is.null(names)) names(param) <- names
      if (!is.null(param)) self$set_par(private$convert_vec2par(param));
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, 2, private$entry_)
      cat("llk:", llk, "\n")
      return(llk)
    },
    
    fit = function(ini_par = NULL, nlm.args = NULL) {
      if (!is.null(ini_par)) self$set_par(ini_par)
      par <- self$par()
      w_par <- private$convert_par2vec(par)
      if (private$print_) cat("Fitting model..........\n")
      t0 <- Sys.time()
      if (is.null(nlm.args)) nlm.args <- list(stepmax = 10)
      args <- c(list(private$calc_negllk, w_par, names = names(w_par), hessian = TRUE), nlm.args)
      mod <- do.call(nlm, args)
      t1 <- Sys.time()
      difft <- t1 - t0
      if (private$print_) cat("Completed model fitting in", difft, attr(difft, "units"), "\n")
      code <- mod$code
      if (code > 2) warning("model failed to converge with nlm code ", code, "\n")
      if (private$print_ & code < 3) cat("Checking convergence.......converged", "\n")
      mle <- mod$estimate
      names(mle) <- names(w_par)
      mle <- private$convert_vec2par(mle)
      self$set_par(mle)
      private$mle_ <- mle
      private$llk_ <- -mod$minimum
      private$V_ <- solve(mod$hessian)
      if (any(diag(private$V) <= 0)) {
        cat("Variance estimates not reliable, do a bootstrap?")
        return(0)
      } else {
        sd <- sqrt(diag(private$V_))
        names(sd) <- names(w_par)
        private$make_results()
      }
      return(invisible())
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
  
  simulate = function(N = NULL, seed = NULL) {
    if (!is.null(N)) N <- self$data()$n()
    new_dat <- simulate_cjs_openscr(self$par(), 
                                    N, 
                                    self$data()$n_occasions(), 
                                    self$data()$traps(), 
                                    self$data()$mesh(), 
                                    move = FALSE, 
                                    time = self$data()$time(), 
                                    seed = seed, 
                                    print = private$print_)
    return(new_dat)
  }, 
  
  par = function() {return(private$par_)},
  mle = function() {return(private$mle_)},
  data = function() {return(private$data_)}, 
  entry = function() {return(private$entry_)}, 
  detectfn = function() {return(private$detfn_)}, 
  
  estimates = function() {
      ests <- NULL
      if (is.null(private$mle_)) {
        ests <- "Fit model using $fit method"
      } else {
        ests$par <- private$results_
      }
      return(ests)
    },
  
    cov_matrix = function() {return(private$V_)}, 
    mle_llk = function() {return(private$llk_)}
),
                   
  private = list(
    data_ = NULL,
    detfn_ = NULL, 
		entry_ = NULL, 
    form_ = NULL, 
    par_ = NULL, 
    link2response_ = NULL, 
    response2link_ = NULL, 
    mle_ = NULL,
    results_ = NULL,
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05, 
		print_ = NULL, 
		detectfn_ = NULL, 
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_det_par <- private$detfn_$npars()
      private$par_ <- vector(mode = "list", length = n_det_par + 1)
      n_par <- numeric(n_det_par)
      for (par in 1:(n_det_par + 1)) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        if (par %in% c(n_det_par + 1) & "t" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("t", private$data_$n_occasions() - 1, sep ="")]
        } else if (par %in% c(n_det_par + 1) & "primary" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("primary", private$data_$n_occasions() - 1, sep ="")]
        }  else {
          n_par[par] <- ncol(X)
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)
        }
        private$par_[[par]] <- par_vec
      }
      names(private$par_) <- names(private$form_)
      return(invisible())
    }, 
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[i]][1] <- do.call(private$response2link_[[i]], 
                                        list(start[[i]]))
      }
      private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
      return(invisible())
    }, 
    
    make_results = function(nsims) {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      par <- private$convert_par2vec(private$mle_) 
      sd <- sqrt(diag(private$V_))
      ci <- private$calc_confint()
      results <- cbind(par, sd, ci$LCL, ci$UCL)
      colnames(results) <- c("Estimate", "Std. Error", "LCL", "UCL")
      rownames(results) <- names(par)
      private$results_ <- results 
      return(invisible())
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
      n_det_par <- self$detectfn()$npars()
      parnames <- self$detectfn()$pars()
      par <- vector(mode = "list", length = n_det_par)
      for (i in 1:n_det_par) {
        par[[i]] <- vec[grep(parnames[i], names)]
        names(par[[i]]) <- gsub(paste0(parnames[i],"."), "", names(par[[i]]))
      }
      names(par) <- parnames 
      par$phi <- vec[grep("phi", names)]
      names(par$phi) <- gsub("phi.", "", names(par$phi))
      return(par)
    }, 
  
   qlog = function(p) {
     q <- p
     q <- ifelse(abs(1 - q) < 1e-16, 1 - 1e-16, q)
     q <- ifelse(abs(q) < 1e-16, 1e-16, q)
     return(qlogis(q))
   }
  )                 
)



