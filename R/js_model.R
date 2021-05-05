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

#' Jolly-Seber model class 
#' 
#' @description Jolly-Seber model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item scr_data: a ScrData object 
#'   \item start: a named list of starting values 
#'   \item print (default = TRUE): if TRUE then helpful output is printed to the screen
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k (if j, k omitted then returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses. Note, the model will simulate 
#'    data using this parameter, but will only present inference based on the maximum likelihood
#'    estimates. 
#'  \item set_mle(mle, V, llk): set maximum likleihood for this model with parameters mle, 
#'  covariance matrix V, and maximum likelihood value llk
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_pr_entry(): computes vector with entry j equal to probability of individual unborn up to occasion j
#'  being born just after occasion j
#'  \item calc_tpms(): returns list of transition probability matrix for each occasion 
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit(): fit the model by obtaining the maximum likelihood estimates. Estimates of
#'        density are obtained from parametric boostrap with nsim resamples. 
#'  \item par(): return current parameter of the model 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#' }
#' 
JsModel <- R6Class("JsModel",
                   inherit = ScrModel,
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, statemod = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      private$data_ <- data
      if (print) cat("Reading formulae.......")
      order <- c("phi", "beta", "D")
      private$read_formula(form, detectfn, statemod, order)
      # add parameters other than detection 
      private$par_type_[private$detfn_$npars() + 1] <- "p1ms"
      private$par_type_[private$detfn_$npars() + 2] <- "pconms"
      private$par_type_[private$detfn_$npars() + 3] <- "m"
      names(private$form_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      # make parameter list 
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"), list("pplink"), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"), list("invpplink"), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      private$read_states() 
      if (print) cat("done\n")
      private$print_ = print 
    },
    
    set_par = function(par) {
      private$par_ <- par
    },
    
    calc_initial_distribution = function() {
      nstates <- private$state_$nstates()
      a0 <- self$get_par("beta", k = 1, m = 1, s = 1:nstates)
      n_mesh <- private$data_$n_meshpts()
      delta <- private$state_$delta() 
      pr0 <- matrix(c(1 - sum(a0*delta), a0*delta, 0), nrow = n_mesh, ncol = nstates + 2, byrow = TRUE)
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:n_mesh) * a 
      for (s in 1:(nstates + 1)) pr0[,s] <- pr0[,s] * D 
      return(pr0)
    },
    
    calc_pr_entry = function() {
      n_occasions <- private$data_$n_occasions()  
      nstates <- self$state()$nstates()
      pr_entry <- matrix(0, nr = n_occasions - 1, nc = nstates)
      for (g in 1:nstates) {
        prod <- 1 - self$get_par("beta", k = 1, m = 1, s = g)
        for (k in 1:(n_occasions - 1)) {
          b <- self$get_par("beta", k = k + 1, m = 1, s = g)
          pr_entry[k, g] <- b / prod
          prod <- prod * (1 - pr_entry[k])
        }
      }
      # prevent numerical error causing prob > 1  
      pr_entry[pr_entry > 1] <- 1 
      return(pr_entry)
    }, 
    
    calc_tpms = function() {
      # compute entry probabilities 
      pr_entry <- self$calc_pr_entry()      
      n_occasions <- private$data_$n_occasions()
      n_primary <- private$data_$n_primary() 
      nstates <- self$state()$nstates()
      delta <- self$state()$delta() 
      tpms <- vector("list", length = n_occasions - 1)
      dt <- diff(private$data_$time())
      ind <- nstates + 1 
      for (k in 1:(n_occasions - 1)) {
        Q <- matrix(0, nr = nstates + 1, nc = nstates + 1)
        Q[-ind, -ind] <- self$state()$trm(k)
        G <- matrix(0, nr = nstates + 2, nc = nstates + 2)
        for (s in 1:nstates) {
          psi <- -log(self$get_par("phi", k = k, m = 1, s = s))
          diag(Q)[s] <- diag(Q)[s] - psi
          Q[s, nstates + 1] <- psi
        }
        G[-1, -1] <- expm(Q * dt[k])
        G[1,-c(1, nstates + 2)] <- pr_entry[k,] * delta 
        G[1, 1] <- 1 - sum(G[1, -1])
        tpms[[k]] <- G
      }
      return(tpms)
    }, 
    
    calc_pr_capture = function() {
      n_occasions <- private$data_$n_occasions("all")
      n_primary <- private$data_$n_primary()
      nstates <- self$state()$nstates()
      kstates <- private$known_states_
      S <- private$data_$n_secondary() 
      if (n_primary == 1) {
        n_primary <- n_occasions
        S <- rep(1, n_occasions)
      }
      enc_rate0 <- self$encrate()
      trap_usage <- usage(private$data_$traps())
      n <- private$data_$n()
      n_meshpts <- private$data_$n_meshpts() 
      n_traps <- private$data_$n_traps()
      capthist <- private$data_$capthist()
      imesh <- private$data_$imesh()
      prob <- C_calc_pr_capture(n, 
                                n_occasions, 
                                n_traps, 
                                n_meshpts, 
                                capthist, 
                                enc_rate0, 
                                trap_usage, 
                                nstates, 
                                1, 
                                1, 
                                kstates, 
                                self$data()$detector_type(),
                                n_primary, 
                                S,
                                rep(0, n), 
                                imesh, 
                                private$data_$capij())
      return(prob)
    },
    
    calc_Dpdet = function() {
      # compute probability of zero capture history 
      n_occasions_all <- private$data_$n_occasions("all")
      n_occasions <- private$data_$n_occasions() 
      nstates <- self$state()$nstates()
      n_primary <- private$data_$n_primary()
      S <- private$data_$n_secondary()
      if (n_primary == 1) {
        n_primary <- n_occasions
        S <- rep(1, n_occasions)
      }
      enc_rate <- self$encrate()
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      j <- 0 
      for (prim in 1:n_primary) { 
        pr_empty[[prim]] <- matrix(1, nr = private$data_$n_meshpts(), nc = nstates + 2)
        pr_empty[[prim]][ , -c(1, nstates + 2)] <- 0 
        for (s in 1:S[prim]) { 
          j <- j + 1
          for (g in 1:nstates) {
            pr_empty[[prim]][, g + 1] <- pr_empty[[prim]][, g + 1] - t(trap_usage[, j]) %*% t(enc_rate[[g]][,,j]) 
          }
        }
        for (g in 1:nstates) pr_empty[[prim]][, g + 1] <- exp(pr_empty[[prim]][, g + 1])
      }
      # average over all life histories 
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      pdet <- C_calc_pdet(n_occasions, pr0, pr_empty, tpms, nstates + 2);
      a <- private$data_$cell_area() 
      D <- self$get_par("D", m = 1:private$data_$n_meshpts()) * a 
      Dpdet <- sum(D) - pdet 
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
      # compute transition probability matrices
      nstates <- self$state()$nstates() + 2
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_Dpdet())
      llk <- llk + self$calc_D_llk()
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
    
    nstates = function() {return(self$state()$nstates() + 2)}, 
    
  print = function() {
    options(scipen = 999)
    if (is.null(private$mle_)) {
      print("Fit model using $fit method")
    } else {
      cat("PARAMETER ESTIMATES (link scale)\n")
      print(signif(private$results_, 4))
    
      cat("--------------------------------------------------------------------------------")
      cat("\n DENSITY (response scale) \n")
      print(signif(private$D_tab_, 4))
      cat("--------------------------------------------------------------------------------")
    }
    options(scipen = 0)
  }
  
),
                   
  private = list(
    Dk_ = NULL, 
    var_ = NULL, 
    confint_ = NULL, 
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                               list(start[[names[i]]]))
      }
      private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                     list(start$phi))
      private$par_$beta[1] <- log(-log(start$beta) / sum(diff(private$data_$time())))
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                list(start$D))
      # compute initial parameters for each jkm
      private$compute_par()
      return(invisible())
    }, 
    
    read_states = function() {
      nstates <- self$state()$nstates() + 2
      kstates <- array(1, dim = c(private$data_$n(), private$data_$n_occasions("all"), nstates))
      covtypes <- private$data_$get_cov_list()$cov_type
      snms <- self$state()$names()
      grpnms <- self$state()$groupnms()
      if ("dead" %in% grpnms) stop("Cannot have a state variable named 'dead'. This is a reserved word.")
      grps <- self$state()$groups()
      alive_states <- 2:(nstates - 1)
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
                if (any(occu)) kstates[i, k, alive_states][!occu] <- -1
              }
            }
          }
          if ("dead" %in% names(s)) {
            if(!is.na(s$dead)) kstates[i, k, -(nstates + 2)] <- -1
          } 
        }
      }
      private$known_states_ <- kstates
      invisible()
    }, 
    
    make_results = function() {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      if (private$print_) cat("Inferring density..........")
      private$infer_D()
      if (private$print_) cat("done\n")
      if (private$print_) cat("Computing variances..........")
      private$calc_var()
      if (private$print_) cat("done\n")
      if (private$print_) cat("Computing confidence intervals..........")
      private$calc_confint()
      if (private$print_) cat("done\n")
      results <- cbind(private$confint_$est$est, private$confint_$est$sds, private$confint_$est$lcl, private$confint_$est$ucl)
      colnames(results) <- c("Estimate", "SE", "LCL", "UCL")
      rownames(results) <- rownames(private$confint_$est) 
      private$results_ <- results 
      
      ## D tab
      D_tab <- matrix(0, nr = self$data()$n_occasions(), nc = 4)
      colnames(D_tab) <- c("Estimate", "SE", "LCL", "UCL")
      time <- private$data_$time() 
      if (private$data_$n_primary() > 1) time <- 1:private$data_$n_primary()
      rownames(D_tab) <- time
      D_tab[, 1] <- private$Dk_
      D_tab[ ,2] <- sqrt((exp(private$var_$Dkvar) - 1) * private$Dk_^2 * exp(private$var_$Dkvar))
      D_tab[, 3] <- private$confint_$Dk$lcl
      D_tab[, 4] <- private$confint_$Dk$ucl

      private$D_tab_ <- D_tab
      
      return(invisible())
    
  }, 
  
  calc_alpha = function(par = NULL, k = NULL) {
    save_par <- NULL
    if (!is.null(par)) {
      save_par <- self$par() 
      self$set_par(private$convert_vec2par(par))
    }
    n_occasions <- self$data()$n_occasions()
    tpms <- self$calc_tpms()
    alpha <- rep(0, n_occasions)
    nstates <- self$state()$nstates()
    delta <- self$state()$delta()
    a0 <- self$get_par("beta", k = 1, m = 1, s = 1:nstates)
    pr <- c(1 - sum(a0*delta), a0*delta, 0)
    alivecols <- 2:(2 + nstates - 1)
    alpha[1] <- sum(pr[alivecols])
    K <- n_occasions
    if (!is.null(k)) K <- k
    if (K > 1) {
      for (occ in 2:K) {
        pr <- pr %*% tpms[[occ - 1]]
        alpha[occ] <- sum(pr[alivecols])
      }
    }
    if (!is.null(k)) alpha <- alpha[K]
    if (!is.null(save_par)) self$set_par(save_par)
    return(alpha)
  }, 
  
  calc_wpdet = function(par = NULL) {
    save_par <- self$par()
    statepar <- self$state()$par() 
    if (!is.null(par)) {
      slen <- length(self$state()$par())
      param2 <- par 
      if (slen > 0) {
        ind <- seq(length(par) - slen + 1, length(par))
        self$state()$set_par(par[ind])
        param2 <- par[-ind]
      }
      self$set_par(private$convert_vec2par(param2));
    }
    pdet <- self$calc_pdet() 
    self$set_par(save_par)
    self$state()$set_par(statepar)
    return(pdet)
  },
  
  calc_var = function(Dk = NULL) {
    n_occasions <- self$data()$n_occasions()
    wpar <- private$convert_par2vec(self$par())
    np <- length(wpar)
    pdet <- self$calc_pdet()
    nparD <- length(wpar[grep("D", names(wpar))])
    exc <- (np-nparD + 1):np
    wpar <- c(wpar, self$state()$par())
    del_pdet <- grad(private$calc_wpdet, wpar)[-exc] 
    # get covariance matrix 
    V <- self$cov_matrix()
    sds <- sqrt(diag(V))
    # theta variance
    dist <- self$data()$distances()
    V_theta <- V[-exc, -exc] * pdet 
    V_theta <- V_theta * self$data()$n()
    # log(D) variance
    Dvar <- t(del_pdet) %*% V_theta %*% del_pdet / pdet^3 + 1 / pdet
    Dvar <- Dvar / (self$data()$area() * self$get_par("D"))
    # if Dk is supplied then compute variance of Dks
    Dk_var <- rep(0, 1)
    Dk <- private$Dk_
    alpha <- private$calc_alpha()
    if (!is.null(Dk)) {
      Dk_var <- numeric(n_occasions) 
      for (k in 1:n_occasions) {
        predictfn <- function(v) {
          slen <- length(self$state()$par())
          param2 <- v 
          if (slen > 0) {
            ind <- seq(length(v) - slen + 1, length(v))
            self$state()$set_par(v[ind])
            param2 <- v[-ind]
          }
          self$set_par(private$convert_vec2par(param2));
          private$infer_D()
          return(log(private$Dk_))
        }
        oldpar <- self$par()
        oldpar2 <- self$state()$par()
        parvec <- private$convert_par2vec(self$par())
        parvec <- c(parvec, self$state()$par())
        J <- jacobian(predictfn, parvec)
        self$set_par(oldpar)
        self$state()$set_par(oldpar2)
        private$infer_D()
        VD <- J %*% V %*% t(J)
        Dk_var <- diag(VD)
      }
    }
    private$var_ <- list(sds = sds, Dvar = Dvar, Dkvar = Dk_var)
    return(invisible())
  }, 
  
   calc_confint = function() {
     V <- private$V_ 
     sds <- sqrt(diag(V))
     est <- private$convert_par2vec(private$mle_)
     slen <- length(self$state()$par()) 
     if (slen > 0) est <- c(est, self$state()$par())
     lev <- 1 - private$sig_level_ / 2 
     alp <- qnorm(lev)
     lcl <- est - alp * sds
     ucl <- est + alp * sds 
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
      estmat <- data.frame(est = est, sds = sds, lcl = lcl, ucl = ucl)
      # Dk, if supplied 
      Dkmat <- NULL
      Dk <- private$Dk_
      if (!is.null(Dk)) {
        Dk_C <- exp(qnorm(lev) * sqrt(private$var_$Dkvar))
        Dk_lcl <- Dk / Dk_C
        Dk_ucl <- Dk * Dk_C
        Dkmat <- data.frame(Dk = Dk, lcl = Dk_lcl, ucl = Dk_ucl) 
      }
      private$confint_ <- list(est = estmat, Dk = Dkmat)
      return(invisible())
    },
  
     infer_D = function() {
       tpms <- self$calc_tpms()
       nstates <- self$state()$nstates()
       delta <- self$state()$delta()
       a0 <- self$get_par("beta", k = 1, m = 1, s = 1:nstates)
       pr0 <- c(1 - sum(a0*delta), a0*delta, 0)
       D <- self$get_par("D")
       alpha <- private$calc_alpha() 
       private$Dk_ <- D * alpha 
       return(invisible())
    }
  )                 
)



