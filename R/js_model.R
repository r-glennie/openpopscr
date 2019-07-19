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
#'  \item simulate(): simulate ScrData object from fitted model
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
      private$form_[[private$detfn_$npars() + 1]] <- form[par_names == "phi"][[1]]
      private$form_[[private$detfn_$npars() + 2]] <- form[par_names == "beta"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"), list("exp"), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"), list("log"), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi", "beta", "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$print_ = print 
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
      if (name == "phi" | name == "beta") j <- 1 
      if (name == "D") {
        j <- 1 
        k <- 1 
      } else {
        m <- 1 
      }
      n_primary <- private$data_$n_primary() 
     if (name == "beta") {
       j <- 1 
       k_save <- k
       k <- 2:(private$data_$n_occasions()) 
       if (n_primary > 1) k <- match(2:n_primary, private$data_$primary())
     } 
     if (name == "phi" & is.null(k)) {
       j <- 1 
       k <- 1:(private$data_$n_occasions() - 1) 
       n_primary <- private$data_$n_primary() 
       if (n_primary > 1) k <- match(2:n_primary, private$data_$primary())
     }
     covs <- private$data_$covs(j = j, k = k, m = m)
     i_par <- which(names(private$form_) == name) 
     if (length(i_par) == 0) stop("No parameter with that name.")
     X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
     if (name %in% c("phi", "beta") & "t" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("t1", sep ="")] 
     }
     if (name %in% c("phi", "beta") & "primary" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("primary1", sep ="")] 
     }
     theta <- private$par_[[i_par]]
     l_par <- which(names(private$link2response_) == name)
     resp <- do.call(private$link2response_[[l_par]], list(X %*% theta))
     if (name == "beta") {
       resp <- c(1, resp)
       # if (!("t" %in% all.vars(private$form_[[i_par]]))) {
       #   dt <- diff(private$data_$time())
       #   resp[-1] <- resp[-1] * dt / sum(dt)
       # } 
       resp <- resp / sum(resp)
       if(!is.null(k_save)) resp <- resp[k_save]
     }
     if (name == "D" & is.null(m)) return(mean(resp))
     return(resp)
    }, 
    
    set_par = function(par) {
      private$par_ <- par
    },
    
    calc_initial_distribution = function() {
      a0 <- self$get_par("beta", j = 1, k = 1)
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1 - a0, a0, 0), nrow = n_mesh, ncol = 3, byrow = TRUE)
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:n_mesh) * a 
      pr0[,1] <- pr0[,1] * D
      pr0[,2] <- pr0[,2] * D
      return(pr0)
    },
    
    calc_pr_entry = function() {
      n_occasions <- private$data_$n_occasions()  
      pr_entry <- rep(0, n_occasions - 1)
      prod <- 1 - self$get_par("beta", j = 1, k = 1, m = 1)
      for (k in 1:(n_occasions - 1)) {
        b <- self$get_par("beta", j = 1, k = k + 1, m = 1)
        pr_entry[k] <- b / prod
        prod <- prod * (1 - pr_entry[k])
      }
      return(pr_entry)
    }, 
    
    calc_tpms = function() {
      # compute entry probabilities 
      pr_entry <- self$calc_pr_entry()      
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
        tpms[[k]] <- matrix(c(1 - pr_entry[k], 0, 0, 
                            pr_entry[k], phi, 0,
                            0, 1 - phi, 1), nrow = 3, ncol = 3)
      }
      return(tpms)
    }, 
    
    calc_pr_capture = function() {
      n_occasions <- private$data_$n_occasions("all")
      n_primary <- private$data_$n_primary()
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
                                3, 
                                self$data()$detector_type(),
                                n_primary, 
                                S,
                                rep(0, n))
      return(prob)
    },
    
    calc_Dpdet = function() {
      # compute probability of zero capture history 
      n_occasions_all <- private$data_$n_occasions("all")
      n_occasions <- private$data_$n_occasions() 
      n_primary <- private$data_$n_primary()
      S <- private$data_$n_secondary()
      if (n_primary == 1) {
        n_primary <- n_occasions
        S <- rep(1, n_occasions)
      }
      enc_rate <- self$calc_encrate()
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      j <- 0 
      for (prim in 1:n_primary) { 
        pr_empty[[prim]] <- matrix(1, nr = private$data_$n_meshpts(), nc = 3)
        pr_empty[[prim]][ , 2] <- 0 
        for (s in 1:S[prim]) { 
          j <- j + 1
          pr_empty[[prim]][, 2] <- pr_empty[[prim]][, 2] - t(trap_usage[, j]) %*% enc_rate[j,,]
        }
        pr_empty[[prim]][,2] <- exp(pr_empty[[prim]][,2])
      }
      # average over all life histories 
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      pdet <- C_calc_pdet(n_occasions, pr0, pr_empty, tpms, 3);
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, 3, rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_Dpdet())
      llk <- llk + self$calc_D_llk()
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
    
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
  }, 
  
  simulate = function(seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    new_dat <- simulate_js_openscr(self$par(), 
                                   self$data()$n_occasions(), 
                                   self$data()$traps(), 
                                   self$data()$mesh(), 
                                   move = FALSE, 
                                   time = self$data()$time(), 
                                   seed = seed, 
                                   print = private$print_)
    return(new_dat)
  }, 
  
  sample_D = function(nsims = 99) {
    if (is.null(private$mle_)) stop("Fit model using $fit method.")
    sds <- sqrt(diag(private$V_))
    return(private$infer_D(nsims, extract_samples = 1))
  }, 
  
  sample_R = function(nsims = 99) {
    if (is.null(private$mle_)) stop("Fit model using $fit method.")
    sds <- sqrt(diag(private$V_))
    return(private$infer_D(nsims, extract_samples = 2))
  } 
  
),
                   
  private = list(
    Dk_ = NULL, 
    var_ = NULL, 
    confint_ = NULL, 
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_det_par <- private$detfn_$npars()
      private$par_ <- vector(mode = "list", length = n_det_par + 1)
      n_par <- numeric(n_det_par + 3)
      for (par in 1:(n_det_par + 3)) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        if (par %in% c(n_det_par+1, n_det_par+2) & "t" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("t", private$data_$n_occasions() - 1, sep ="")]
        } else if (par %in% c(n_det_par+1, n_det_par+2) & "primary" %in% all.vars(private$form_[[par]])) {
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
      private$par_$beta[1] <- do.call(private$response2link_$beta,
                                      list(c(1 / start$beta - 1)/self$data()$n_occasions()))
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                list(start$D / private$data_$n_meshpts()))
      return(invisible())
    }, 
    
    make_results = function() {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      par <- private$convert_par2vec(private$mle_) 
      if (private$print_) cat("Inferring density..........")
      private$infer_D()
      if (private$print_) cat("done\n")
      if (private$print_) cat("Computing variances..........")
      private$calc_var()
      if (private$print_) cat("done\n")
      if (private$print_) cat("Computing confidence intervals..........")
      private$calc_confint()
      if (private$print_) cat("done\n")
      results <- cbind(par, private$var_$sds, private$confint_$est$lcl, private$confint_$est$ucl)
      colnames(results) <- c("Estimate", "SE", "LCL", "UCL")
      rownames(results) <- names(par)
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
    a0 <- self$get_par("beta", k = 1)
    pr <- c(1 - a0, a0, 0)
    alpha[1] <- pr[2]
    K <- n_occasions
    if (!is.null(k)) K <- k
    if (K > 1) {
      for (occ in 2:K) {
        pr <- pr %*% tpms[[k - 1]]
        alpha[occ] <- pr[2]
      }
    }
    if (!is.null(k)) alpha <- alpha[K]
    if (!is.null(save_par)) self$set_par(save_par)
    return(alpha)
  }, 
  
  calc_wpdet = function(par = NULL) {
    save_par <- self$par() 
    if (!is.null(par)) self$set_par(private$convert_vec2par(par))
    pdet <- self$calc_pdet() 
    self$set_par(save_par)
    return(pdet)
  },
  
  calc_var = function(Dk = NULL) {
    n_occasions <- self$data()$n_occasions()
    wpar <- private$convert_par2vec(self$par())
    np <- length(wpar)
    pdet <- self$calc_pdet()
    nparD <- length(wpar[grep("D", names(wpar))])
    exc <- (np-nparD + 1):np
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
    if (!is.null(Dk)) {
      Dk_var <- numeric(n_occasions) 
      for (k in 1:n_occasions) {
        del_alpha <- grad(private$calc_alpha, wpar, k = k)[-exc]
        sig_alpha <- t(del_alpha) %*% V_theta %*% del_alpha
        Dk_var[k] <- Dk[k] * Dvar / self$get_par("D") + sig_alpha / self$get_par("D")
      }
      Dk_var <- Dk_var / (self$data()$area() * Dk)
    }
    private$var_ <- list(sds = sds, Dvar = Dvar, Dkvar = Dk_var)
    return(invisible())
  }, 
  
   calc_confint = function() {
     V <- private$V_ 
     sds <- sqrt(diag(V))
     est <- private$convert_par2vec(private$mle_)
     lev <- 1 - private$sig_level_ / 2 
     alp <- qnorm(lev)
     lcl <- est - alp * sds
     ucl <- est + alp * sds 
      estmat <- data.frame(lcl = lcl, ucl = ucl)
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
       a0 <- self$get_par("beta", k = 1)
       pr0 <- c(1 - a0, a0, 0)
       private$Dk_ <- C_calc_D(self$get_par("D"), self$data()$n_occasions(), pr0, tpms)
       return(invisible())
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
      par$beta <- vec[grep("beta", names)]
      names(par$beta) <- gsub("beta.", "", names(par$beta))
      par$D <- vec[grep("D", names)]
      names(par$D) <- gsub("D.", "", names(par$D))
      return(par)
    }
  )                 
)



