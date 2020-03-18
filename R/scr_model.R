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
#'  \item predict(name, newdata, type, se, alpha): returns predictions for named parameter 
#'        newdata can be an alternative design matrix for that parameter, 
#'        type = "response" provides predictions on response scale, otherwise link scale 
#'        se = TRUE leads to computation of confidence intervals based on delta method 
#'        alpha = 0.95 is the confidence level 
#'        nsims = number of simulations used for mlogit confidnce intervals 
#'        seed = seed used for simulations 
#'  \item pr_state: return a list (one entry for each individual) of arrays with entry (m, s, k) being
#'        the probability the individual was at mesh point m, in state s, in occasion k given the observed
#'        data and current parameters
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
#'  \item fit(): fit the model by obtaining the maximum likelihood estimates
#'  \item set_mle(par, V, llk): set par to maximum likelihood par with variance matrix V and value llk
#'  \item par(): return current parameter of the model 
#'  \item npar(): number of parameters in the model 
#'  \item state(): returns stored StateModel object 
#'  \item design_mat(): returns a list of the design matrices used for each parameter 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to
#'  \item set_par(): overwrite parameter stored, current value returned by $par() 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#' }
#' 
ScrModel <- R6Class("ScrModel", 
  public = list(
    
    initialize = function(form, data, start, kbuf = 0, detectfn = NULL, statemod = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      private$data_ <- data
      if (print) cat("Reading formulae.......")
      private$read_formula(form, detectfn, statemod)
      # add parameters other than detection 
      private$par_type_[private$detfn_$npars() + 1] <- "m"
      names(private$form_) <- c(private$detfn_$pars(), "D")
      # make parameter list
      private$make_par() 
      # set link functions and inverse link functions
      private$link2response_ <- c(private$detfn_$link2response(), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      private$read_states() 
      if (print) cat("done\n")
      private$print_ = print 
      private$kbuf_ = kbuf
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL, s = NULL) {
      if (!(name %in% names(private$computed_par_))) stop("No parameter with that name.")
      ipar <- match(name, names(private$computed_par_))
      type <- private$par_type_[ipar]
      if (is.null(j)) j <- 1:private$data_$n_traps() 
      if (is.null(k)) {
        if (type == "k1ms" | type == "k1m") {
          k <- 1:(private$data_$n_occasions("all") - ifelse(private$data_$n_primary() > 1, private$data_$n_secondary()[private$data_$n_primary()], 1))
        } else if (type == "p1ms" | type == "p1m") {
          k <- 1:(private$data_$n_occasions() - 1) 
        } else if (type == "pconms") {
          k <- 1:private$data_$n_occasions()
        } else {
          k <- 1:private$data_$n_occasions("all") 
        }
      }
      if (is.null(s)) s <- 1
      if (!identical(private$par_, private$last_par_)) private$compute_par() 
      if (type == "jk") {
        return(private$computed_par_[[ipar]][k, j])
      } else if (type == "jks") {
        return(private$computed_par_[[ipar]][k, j, s])
      } else if (type == "km" | type == "k1m" | type == "kconm") {
        j <- 1
        mnew <- m 
        if (is.null(mnew)) mnew <- 1:private$data_$n_meshpts()
        res <- private$computed_par_[[ipar]][k, mnew]
        if (is.null(m) & length(k) > 1) return(rowMeans(res))
        if (is.null(m) & length(k) == 1) return(mean(res))
        return(res)
      } else if (type == "k1ms" | type == "kconms") {
        mnew <- m 
        if (is.null(mnew)) mnew <- 1:private$data_$n_meshpts()
        res <- private$computed_par_[[ipar]][k, mnew, s]
        return(res)
      } else if (type == "p1ms" | type == "pconms") {
          mnew <- m 
          if (is.null(mnew)) mnew <- 1:private$data_$n_meshpts()
          res <- private$computed_par_[[ipar]][k, mnew, s]
          return(res)
      } else if (type == "m") {
        j <- 1 
        k <- 1
        mnew <- m 
        if (is.null(mnew)) mnew <- 1:private$data_$n_meshpts() 
        res <- private$computed_par_[[ipar]][mnew]
        if (is.null(m)) return(mean(res))
        return(res)
      }
    }, 
    
    predict = function(name, newdata = NULL, type = "response", se = FALSE, alpha = 0.95, nsims = 100, seed = 15185) {
      set.seed(seed)
      on.exit(set.seed(round(runif(1, 0, 10000))))
      if (is.null(private$mle_)) print("Fit model using $fit method")
      if (!(name %in% names(private$computed_par_))) stop("No parameter with that name.")
      ipar <- match(name, names(private$computed_par_))
      X <- newdata
      if (is.null(X)) {
        X <- private$X_[[ipar]]
      } else {
        if (!is.data.frame(X)) stop("newdata must be a data frame.")
        if (is.data.frame(X)) X <- as.matrix(X)
      }
      if (ncol(X) != ncol(private$X_[[ipar]])) stop("newdata has wrong number of columns.")
      if (!identical(colnames(X)[-1], colnames(private$X_[[ipar]])[-1])) stop("newdata columns in wrong order or wrong names.")
      
      old_par <- private$convert_par2vec(self$mle())
      predicterfn <- function(v) {
        self$set_par(private$convert_vec2par(v));
        return(X %*% private$par_[[ipar]])
      }
      est <- predicterfn(old_par) 
      if (se) {
        g <- jacobian(predicterfn, old_par)
        V <- g %*% private$V_ %*% t(g) 
        sds <- sqrt(diag(V))
        z <- qnorm(1 - (1 - alpha) / 2)
        lcl <- est - z * sds
        ucl <- est + z * sds
      }
      if (type == "response") {
        if (se) {
          if (private$link2response_[[ipar]] != "mlogit") {
            lcl <- as.vector(do.call(private$link2response_[[ipar]], list(lcl)))
            ucl <- as.vector(do.call(private$link2response_[[ipar]], list(ucl)))
          } else {
            sims <- t(rmvn(nsims, est[,1], V))
            sims <- do.call(private$link2response_[[ipar]], list(sims))
            lcl <- apply(sims, 1, FUN = function(x) {quantile(x, prob = (1 - alpha)/2)})
            ucl <- apply(sims, 1, FUN = function(x) {quantile(x, prob = 1 - (1 - alpha)/2)})
          }
        }
        est <- as.vector(do.call(private$link2response_[[ipar]], list(est)))
      }
      res <- est
      if(se) res <- data.frame(est = est, lcl = lcl, ucl = ucl)
      if(se) {attributes(res)$vcov <- V;  attributes(res)$g <- g}
      self$set_par(self$mle());
      return(res)
    }, 
    
    pr_state = function() {
      fw <- private$calc_forwback() 
      nocc <- private$data_$n_occasions("all")
      n <- private$data_$n()
      nmesh <- private$data_$n_meshpts()
      nstates <- self$nstates()
      pred <- vector(mode = "list", length = n)
      for (i in 1:n) {
        pred[[i]] <- array(0, dim = c(nmesh, nstates, nocc))
        for (j in 1:nocc) {
          pred[[i]][,,j] <- exp(fw$lalpha[[i]][,,j] + fw$lbeta[[i]][,,j])
          pred[[i]][,,j] <- pred[[i]][,,j] / sum(pred[[i]][,,j])
        }
      }
      return(pred)
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
      nstates <- private$state_$nstates()
      delta <- private$state_$delta() 
      pr0 <- matrix(delta, nrow = n_mesh, ncol = nstates, byrow = TRUE)
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:n_mesh) * a 
      pr0 <- pr0 * D
      return(pr0)
    },
    
    calc_encrate = function(transpose = FALSE) {
      dist <- private$data_$distances()
      n_occasions <- private$data_$n_occasions("all")
      nstates <- self$state()$nstates()
      enc_rate <- vector(mode = "list", length = nstates)
      n_det_par <- self$detectfn()$npars()
      det_par <- vector(mode = "list", length = n_det_par)
      for (s in 1:nstates) {
        enc_rate[[s]] <- array(0, dim = c(n_occasions, nrow(dist), ncol(dist))) 
        for (k in 1:n_occasions) {
          for (dpar in 1:n_det_par) det_par[[dpar]] <- as.vector(self$get_par(self$detectfn()$par(dpar), k = k, m = 1, s = s))
          if (any(sapply(det_par, anyNA))) {
            if (all(sapply(det_par, anyNA))) {
              enc_rate[[s]][k,,] <- 0
            } else {
              stop("Make sure that either all detection parameters have a state variable with NA's in the same places or none do.")
            }
          } else {
            enc_rate[[s]][k,,] <- self$detectfn()$h(dist, det_par)
          }
        }
        # add epsilon to stop log(0.0)
        enc_rate[[s]] <- enc_rate[[s]] + 1e-16
      }
      if (transpose) {
        for (s in 1:nstates) {
          temp <- enc_rate[[s]]
          enc_rate[[s]] <- array(0, dim = c(ncol(dist), nrow(dist), n_occasions)) 
          for (k in 1:n_occasions) enc_rate[[s]][,,k] <- t(temp[k,,])
        }
      }
     return(enc_rate) 
    }, 
    
    calc_pr_capture = function() {
      n_occasions <- private$data_$n_occasions()
      enc_rate0 <- self$calc_encrate(transpose = TRUE)
      trap_usage <- usage(private$data_$traps())
      n <- private$data_$n()
      n_meshpts <- private$data_$n_meshpts() 
      n_states <- private$state_$nstates() 
      n_traps <- private$data_$n_traps()
      capthist <- private$data_$capthist()
      kstates <- private$known_states_
      dist <- private$data_$distances()
      imesh <- private$data_$imesh()
      prob <- C_calc_pr_capture(n, 
                                n_occasions, 
                                n_traps, 
                                n_meshpts,
                                capthist, 
                                enc_rate0, 
                                trap_usage, 
                                n_states, 
                                0, 
                                0,
                                kstates,
                                self$data()$detector_type(), 
                                n_occasions, 
                                rep(1, n_occasions),
                                rep(0, n), 
                                imesh,
                                private$data_$capij())
      return(prob)
    },
    
    calc_Dpdet = function() {
      enc_rate <- self$calc_encrate() 
      nstates <- self$state()$nstates()
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      for (j in 1:private$data_$n_occasions()) {
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = nstates)
        for (g in 1:nstates) {
          pr_empty[[j]][, g] <- exp(-t(trap_usage[, j]) %*% enc_rate[[g]][j,,])
        }
      }
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      Dpdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms, nstates);
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
      # initial distribution 
      pr0 <- self$calc_initial_distribution()
      # compute probability of capture histories 
      # across all individuals, occasions and traps 
      pr_capture <- self$calc_pr_capture()
      # compute likelihood for each individual
      n <- private$data_$n()
      n_occasions <- private$data_$n_occasions()
      n_meshpts <- private$data_$n_meshpts() 
      # get tpms for state model 
      nstates <- self$state()$nstates() 
      tpms <- self$calc_tpms()
      # compute log-likelihood
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, private$data_$n()))
      llk <- llk - n * log(self$calc_Dpdet())
      llk <- llk + self$calc_D_llk()
      if (private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
    
    fit = function(nlm.args = NULL) {
      par <- self$par()
      # scrmodel parameters 
      w_par <- private$convert_par2vec(par)
      # add in state model parameters
      w_par <- c(w_par, self$state()$par())
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
      if (private$print_) cat("Computing variances.......")
      mle <- mod$estimate
      names(mle) <- names(w_par)
      llk <- -mod$minimum
      V <- solve(mod$hessian)
      if (private$print_) cat("done\n")
      self$set_mle(mle, V, llk)
      return(invisible())
    }, 
    
    set_mle = function(par, V, llk) {
      slen <- length(self$state()$par())
      param2 <- par 
      if (slen > 0) {
        ind <- seq(length(par) - slen + 1, length(par))
        self$state()$set_par(par[ind])
        param2 <- par[-ind]
      }
      mle <- private$convert_vec2par(param2);
      self$set_par(mle)
      private$mle_ <- mle 
      private$llk_ <- llk 
      private$V_ <- V
      if (any(diag(V) <= 0)) warning("Variance estimates not reliable, do a bootstrap.")
      sd <- sqrt(diag(private$V_))
      names(sd) <- names(par)
      private$make_results()
    },
    
    print = function() {
      options(scipen = 999)
      if (is.null(private$mle_)) {
        print("Fit model using $fit method")
      } else {
        cat("PARAMETER ESTIMATES (link scale)\n")
        print(signif(private$results_, 4))
        slen <- length(self$state()$par())
        if (slen > 0) {
          cat("\n")
          self$state()$print()
        }
      }
      options(scipen = 0)
    }, 
    
    calc_tpms = function() {
      noccasions <- private$data_$n_occasions()
      if (noccasions == 1) return(diag(1))
      tpms <- vector("list", length = noccasions - 1)
      dt <- diff(private$data_$time())
      for (k in 1:(noccasions - 1)) {
        tpms[[k]] <- self$state()$tpm(k = k, dt = dt[k])
      }
      return(tpms)
    }, 
    
    par = function() {return(private$par_)},
    npar = function() {return(length(unlist(self$par())) + length(self$state()$par()))}, 
    state = function() {return(private$state_)}, 
    design_mats = function() {return(private$X_)}, 
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
    mle_llk = function() {return(private$llk_)}, 
    nstates = function() {return(self$state()$nstates())}
    
  ),
  
  private = list(
    data_ = NULL,
    detfn_ = NULL, 
    state_ = NULL, 
    known_states_ = NULL, 
    form_ = NULL, 
    par_ = NULL, 
    last_par_ = NULL,
    computed_par_ = NULL, 
    par_type_ = NULL, 
    X_ = NULL, 
    link2response_ = NULL, 
    response2link_ = NULL, 
    mle_ = NULL,
    results_ = NULL,
    D_tab_ = NULL, 
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05, 
    print_  = NULL,
    kbuf_ = NULL, 
    
    read_formula = function(form, detectfn, statemod, order = NULL) {
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$par_type_ <- rep("", length(par_names))
      # state model 
      if (is.null(statemod)) {
        private$state_ <- StateModel$new(private$data_, 
                                         c("Alive"), 
                                         matrix(".", nr = 1, nc = 1), 
                                         list(delta = c(1), trm = matrix(0, nr = 1, nc = 1)))
      } else {
        private$state_ <- statemod$clone()
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
        private$par_type_[i] <- "jks"
      }
      if (!is.null(order)) {
        ord <- c(private$detfn_$pars(), order)
        for (i in 1:length(private$form_)) private$form_[[i]] <- form[par_names == ord[i]][[1]]
      }
      return(invisible())
    }, 
    
    make_par = function() {
      samp_cov_jk <- private$data_$covs(m = 1)
      samp_cov_km <- private$data_$covs(j = 1)
      samp_cov_pm <- private$data_$covs(j = 1)
      samp_cov_m <- private$data_$covs(j = 1, k = 1)
      npar <- length(private$par_type_)
      private$par_ <- vector(mode = "list", length = npar)
      trapocc <- expand.grid(list(occ = 1:private$data_$n_occasions("all"), 
                             trap = 1:private$data_$n_traps()))
      occm <- expand.grid(list(occ = 1:private$data_$n_occasions("all"), 
                               mesh = 1:private$data_$n_meshpts()))
      primm <- expand.grid(list(occ = 1:private$data_$n_occasions(), 
                               mesh = 1:private$data_$n_meshpts()))
      covslst <- private$data_$get_cov_list()
      covs <- covslst$cov_type
      tempdatjk <- tempdatkm <- tempdatm <- tempdatpm <- NULL
      covnms <- names(covslst$cov)
      for (i in 1:length(covs)) {
        if (covs[i] %in% c("i", "ik")) next 
        k <- switch(covs[i], "k" = 1, "j" = 2, "m" = 3, 
                    "kj" = 4, "km" = 5, "p" = 6, "pm" = 7)
        if (k == 1) {
          tempdatjk[[covnms[i]]] <- samp_cov_jk[[i]][trapocc[,k]]
          tempdatkm[[covnms[i]]] <- samp_cov_km[[i]][occm[,k]]
        } else if (k == 2) {
          tempdatjk[[covnms[i]]] <- samp_cov_jk[[i]][trapocc[,k]]
        } else if (k == 3) {
          tempdatkm[[covnms[i]]] <- samp_cov_km[[i]][occm[,k - 1]]
          tempdatm[[covnms[i]]] <- as.vector(samp_cov_m[[i]])
        } else if (k == 4) {
          tempdatjk[[covnms[i]]] <- as.vector(samp_cov_jk[[i]])
        } else if (k == 5) {
          tempdatkm[[covnms[i]]] <- as.vector(samp_cov_km[[i]])
        } else if (k == 6) {
          tempdatpm[[covnms[i]]] <- samp_cov_pm[[i]][primm[,1]]
        } else if (k == 7) {
          tempdatpm[[covnms[i]]] <- as.vector(samp_cov_pm[[i]])
        }
      }
      tempdatjk <- as.data.frame(tempdatjk)
      tempdatkm <- as.data.frame(tempdatkm)
      tempdatpm <- as.data.frame(tempdatpm)
      tempdatm <- as.data.frame(tempdatm)
      if ("par" %in% c(names(tempdatjk), names(tempdatkm), names(tempdatm), names(tempdatpm))) stop("Cannot call covariates 'par'. Change name.")
      tempdatjk$par <- 1 
      tempdatkm$par <- 1
      tempdatpm$par <- 1
      tempdatm$par <- 1 
      # remove covariates from last primary/secondary occasion for k1m parameters
      rm <- seq(private$data_$n_occasions("all"), nrow(tempdatkm), by = private$data_$n_occasions("all"))
      tempdatk1m <- droplevels(tempdatkm[-rm,])
      rm <- seq(private$data_$n_occasions(), nrow(tempdatpm), by = private$data_$n_occasions())
      tempdatp1m <- droplevels(tempdatpm[-rm,])
      # create tempdat for kjs parameters 
      nstates <- self$state()$nstates()
      tempdatjks <- do.call("rbind", replicate(nstates, tempdatjk, simplify = FALSE))
      grps <- self$state()$groups()
      ngrps <- ncol(grps)
      grpnm <- self$state()$groupnms()
      for (g in 1:ngrps) {
        tempdatjks[[grpnm[g]]] <- rep(1:nstates, each = nrow(tempdatjk))
        tempdatjks[[grpnm[g]]] <- grps[tempdatjks[[grpnm[g]]], g]
      }
      # create tempdat for k1ms parameters 
      tempdatk1ms <- do.call("rbind", replicate(nstates, tempdatk1m, simplify = FALSE))
      for (g in 1:ngrps) {
        tempdatk1ms[[grpnm[g]]] <- rep(1:nstates, each = nrow(tempdatk1m))
        tempdatk1ms[[grpnm[g]]] <- grps[tempdatk1ms[[grpnm[g]]], g]
      }
      # create tempdat for p1ms parameters 
      tempdatp1ms <- do.call("rbind", replicate(nstates, tempdatp1m, simplify = FALSE))
      for (g in 1:ngrps) {
        tempdatp1ms[[grpnm[g]]] <- rep(1:nstates, each = nrow(tempdatp1m))
        tempdatp1ms[[grpnm[g]]] <- grps[tempdatp1ms[[grpnm[g]]], g]
      }
      # collate tempdats together 
      tempdat <- list(tempdatjk, tempdatkm, tempdatm, tempdatk1m, tempdatjks, tempdatk1ms, tempdatpm, tempdatp1ms)
      parloc <- sapply(tempdat, FUN = function(x) {min(which(names(x) == "par"))})
      for (par in 1:npar) {
        k <- switch(private$par_type_[par], "jk" = 1, "km" = 2, "m" = 3, "k1m" = 4, "kconm" = 4, "jks" = 5, "k1ms" = 6, "kconms" = 6, "pm" = 7, "p1ms" = 8, "pconms" = 8)
        parname <- as.character(private$form_[[par]][[2]])
        names(tempdat[[k]])[parloc[k]] <- parname
        private$X_[[par]] <- openpopscrgam(private$form_[[par]], data = tempdat[[k]])
        names(private$X_)[[par]] <- parname
        par_vec <- rep(0, ncol(private$X_[[par]]))
        names(par_vec) <- colnames(private$X_[[par]])
        private$par_[[par]] <- par_vec
      }
      names(private$par_) <- names(private$form_)
      return(invisible())
    }, 
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                        list(start[[names[i]]]))
      }
      names(private$par_) <- c(names, "D")
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                list(start$D))
      # compute initial parameters for each jkm
      private$compute_par()
      return(invisible())
    }, 
    
    compute_par = function() {
      private$last_par_ <- private$par_
      npar <- length(private$par_type_)
      private$computed_par_ <- vector(mode = "list", length = npar)
      for (par in 1:npar) {
        dimsize <- switch(private$par_type_[par], "jk" = 2,
                          "jks" = 3, 
                          "km" = 2, 
                          "pm" = 2, 
                          "m" = 2,
                          "k1m" = 2, 
                          "kconm" = 2, 
                          "k1ms" = 3, 
                          "kconms" = 3,
                          "p1ms" = 3, 
                          "pconms" = 3
                          )
        dim <- rep(0, dimsize)
        dim[1] <- switch(private$par_type_[par], "jk" = private$data_$n_occasions("all"),
                         "jks" = private$data_$n_occasions("all"), 
                     "km" = private$data_$n_occasions("all"), 
                     "k1ms" = private$data_$n_occasions("all") - 
                       ifelse(private$data_$n_primary() > 1, private$data_$n_secondary()[private$data_$n_primary()], 1), 
                     "m" = private$data_$n_meshpts(),
                      "k1m" = private$data_$n_occasions("all") - 
                       ifelse(private$data_$n_primary() > 1, private$data_$n_secondary()[private$data_$n_primary()], 1), 
                     "kconm" = private$data_$n_occasions("all") - 
                       ifelse(private$data_$n_primary() > 1, private$data_$n_secondary()[private$data_$n_primary()], 1), 
                     "kconms" = private$data_$n_occasions("all") - 
                       ifelse(private$data_$n_primary() > 1, private$data_$n_secondary()[private$data_$n_primary()], 1), 
                     "pm" = private$data_$n_occasions(), 
                     "p1ms" = private$data_$n_occasions() - 1, 
                     "pconms" = private$data_$n_occasions() - 1)
        if (dimsize > 1) {
        dim[2] <- switch(private$par_type_[par], "jk" = private$data_$n_traps(),
                         "jks" = private$data_$n_traps(), 
                     "km" = private$data_$n_meshpts(), 
                     "k1ms" = private$data_$n_meshpts(), 
                     "m" = 1,
                     "k1m" = private$data_$n_meshpts(),
                      "kconm" = private$data_$n_meshpts(), 
                     "kconms" = private$data_$n_meshpts(), 
                     "pm" =  private$data_$n_meshpts(), 
                     "p1ms" = private$data_$n_meshpts(), 
                     "pconms" = private$data_$n_meshpts())
        } 
        if (dimsize > 2) {
          dim[3] <- switch(private$par_type_[par], "jks" = private$state_$nstates(), "k1ms" = private$state_$nstates(), 
                           "kconms" = private$state_$nstates(), "p1ms" = private$state_$nstates(), 
                           "pconms" = private$state_$nstates())
        }
        private$computed_par_[[par]] <- do.call(private$link2response_[[par]],
                                                list(array(private$X_[[par]] %*% private$par_[[par]], 
                                                     dim = dim)))
      }
      names(private$computed_par_) <- names(private$form_)
      return(invisible())
    }, 
    
    read_states = function() {
      kstates <- array(1, dim = c(private$data_$n(), private$data_$n_occasions("all"), self$state()$nstates()))
      covtypes <- private$data_$get_cov_list()$cov_type
      snms <- self$state()$names()
      grpnms <- self$state()$groupnms()
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
        }
      }
      private$known_states_ <- kstates
      invisible()
    }, 
    
    make_results = function() {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      par <- private$convert_par2vec(private$mle_)
      sd <- sqrt(diag(private$V_))
      if (private$print_) cat("Computing confidence intervals.......")
      ci <- private$calc_confint()
      if (private$print_) cat("done\n")
      results <- cbind(ci$est, ci$sds, ci$LCL, ci$UCL)
      colnames(results) <- c("Estimate", "Std. Error", "LCL", "UCL")
      rownames(results) <- names(par)
      private$results_ <- results 
      return(invisible())
      
    }, 
    
    calc_confint = function() {
      V <- private$V_ 
      sds <- sqrt(diag(V))
      est <- private$convert_par2vec(private$mle_)
      slen <- length(self$state()$par()) 
      if (slen > 0) est <- c(est, self$state()$par())
      alp <- qnorm(1 - private$sig_level_ / 2)
      lcl <- est - alp * sds
      ucl <- est + alp * sds 
      names(lcl) <- names(ucl) <- names(est)
      if (slen > 0) {
        ind <- seq(length(est) - slen + 1, length(est))
        slcl <- lcl[ind]
        sucl <- ucl[ind]
        sest <- est[ind]
        ssds <- sds[ind]
        self$state()$calc_confint(sest, ssds, slcl, sucl)
        est <- est[-ind]
        sds <- sds[-ind]
        lcl <- lcl[-ind]
        ucl <- ucl[-ind]
      }
      return(list(est = est, sds = sds, LCL = lcl, UCL = ucl))
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
      nstates <- self$nstates()
      tpms <- self$calc_tpms()
      # compute forward-backward 
      if (forw) lalpha <- C_calc_alpha(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, private$data_$n()))
      if (back) lbeta <- C_calc_beta(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, private$data_$n()))
      if (forw & back) return(list(lalpha = lalpha, lbeta = lbeta))
      if (forw) return(lalpha)
      if (back) return(lbeta)
      return(0)
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
      names <- names(vec)
      npar <- length(private$par_)
      parnames <- names(private$par_)
      par <- vector(mode = "list", length = npar)
      for (i in 1:npar) {
        par[[i]] <- vec[grep(parnames[i], names)]
        names(par[[i]]) <- gsub(paste0(parnames[i],"."), "", names(par[[i]]))
      }
      names(par) <- parnames 
      return(par)
    }, 
    
    check_input = function(form, data, start, detectfn, print) {
      if (!is.list(form) | any(sapply(form, FUN = class) != "formula")) stop("Invalid list of formulae.")
      if (!("ScrData" %in% class(data))) stop("Must supply a ScrData object for data.")
      if (!is.list(start)) stop("start must be a list of starting values for each parameter.")
      if (!is.null(detectfn)) {
        if (!("DetectFn" %in% class(detectfn))) stop("detectfn must be a DetectFn object.")
      } 
      if (!is.logical(print)) stop("print must be TRUE or FALSE.")
      return(0)
    }
  )
)


