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
#'   \item detectfn: either a character of "HHN" (hazard half-normal, default), "HN" (half-normal),
#'         or an object of class DetectFn
#'   \item print (default = TRUE): if TRUE then helpful output is printed
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k at mesh point m, similar to ScrData$covs() function
#'  \item calc_D_llk(): computes the likelihood contribution from the point process model for D
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead) and mesh
#'  \item calc_encrate(transpose = FALSE): compute encounter rate for each mesh x trap x occasion, can be 
#'        return transposed as occasion x mesh x trap
#'  \item calc_pr_capture(): returns list of arrays, one per individual, with (m, s, k) entry being 
#'        the probability of the capture record on occasion k given activity centre is at mesh point m and 
#'        individual is in life state s
#'  \item calc_Dpdet(): compute the integral int D(x)p(x) dx <- the overall intensity of detected individuals in the 
#'        survey area 
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit: fit the model by obtaining the maximum likelihood estimates 
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
    
    initialize = function(form, data, start, detectfn = NULL, print = TRUE) {
      private$data_ <- data
      if (print) cat("Reading formulae.......")
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      if (!("D" %in% par_names)) {
        private$form_ <- c(private$form_, list(D ~ 1))
        par_names <- c(par_names, "D") 
      }
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
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c(private$detfn_$pars(), "D")
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$print_ = print 
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
      if (name == "D") {
        j <- 1 
        k <- 1 
      } else {
        m <- 1 
      }
      covs <- private$data_$covs(j = j, k = k, m = m)
      i_par <- which(names(private$form_) == name) 
      if (length(i_par) == 0) stop("No parameters with that name.")
      X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
      theta <- private$par_[[i_par]]
      l_par <- which(names(private$link2response_) == name)
      resp <- do.call(private$link2response_[[l_par]], list(X %*% theta))
      if (name == "D" & is.null(m)) return(mean(resp))
      return(resp)
    }, 
    
    calc_D_llk = function() {
      n <- private$data_$n()
      Dpdet <- self$calc_Dpdet()
      llk <- n * log(sum(Dpdet)) - sum(Dpdet) - lfactorial(n)
      names(llk) <- NULL 
      return(llk)
    },
    
    calc_initial_distribution = function() {
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(1, nrow = n_mesh, ncol = 1, byrow = TRUE)
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:n_mesh) * a 
      pr0 <- pr0 * D
      return(pr0)
    },
    
    calc_encrate = function(transpose = FALSE) {
      dist <- private$data_$distances()
      if (transpose) dist <- t(dist)
      n_occasions <- private$data_$n_occasions("all")
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
      n_occasions <- private$data_$n_occasions()
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
                                1, 
                                self$data()$detector_type(), 
                                n_occasions, 
                                rep(1, n_occasions),
                                rep(0, n))
      return(prob)
    },
    
    calc_Dpdet = function() {
      enc_rate <- self$calc_encrate() 
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      for (j in 1:private$data_$n_occasions()) {
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = 1)
        pr_empty[[j]][, 1] <- exp(-t(trap_usage[, j]) %*% enc_rate[j,,])
      }
      pr0 <- self$calc_initial_distribution()
      tpms <- list(matrix(0, nr = 2, nc = 2))
      Dpdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms, 1);
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:private$data_$n_meshpts()) * a 
      Dpdet <- sum(D) - Dpdet
      return(Dpdet)
    },
    
    calc_pdet = function() {
      savepar <- self$par()
      newpar <- self$par() 
      newpar$D <- rep(0, length(savepar$D))
      newpar$D[1] <- log(1.0 / self$data()$area())
      self$set_par(newpar)
      pdet <- self$calc_Dpdet()  
      self$set_par(savepar)
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, 1, rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_Dpdet())
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
      if (is.null(nlm.args)) nlm.args <- list(stepmax = 10)
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
    
    par = function() {return(private$par_)},
    mle = function() {return(private$mle_)},
    data = function() {return(private$data_)}, 
    detectfn = function() {return(private$detfn_)}, 
    
    set_par = function(par) {
      private$par_ <- par
    },
    
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
    detfn_ = NULL, 
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
    print_  = NULL, 
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_det_par <- private$detfn_$npars()
      private$par_ <- vector(mode = "list", length = n_det_par + 1)
      n_par <- numeric(n_det_par + 1)
      for (par in 1:(n_det_par + 1)) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        par_vec <- rep(0, n_par[par])
        names(par_vec) <- colnames(X)
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
      names(private$par_) <- c(names, "D")
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                list(start$D / private$data_$n_meshpts()))
      return(invisible())
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
      parnames <- c(self$detectfn()$pars(), "D")
      par <- vector(mode = "list", length = n_det_par)
      for (i in 1:(n_det_par + 1)) {
        par[[i]] <- vec[grep(parnames[i], names)]
        names(par[[i]]) <- gsub(paste0(parnames[i],"."), "", names(par[[i]]))
      }
      names(par) <- parnames 
      return(par)
    }
  )
)


