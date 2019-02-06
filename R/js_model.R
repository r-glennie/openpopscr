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
#'   \item num_cores (optional, default = 1): number of processors cores to use 
#'   in parallelised code 
#'   \item print (defualt = TRUE): if TRUE then helpful output is printed to the screen
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
  public = list(
    
    initialize = function(form, data, start, num_cores = 1, print = TRUE) {
      private$data_ <- data
      if (print) cat("Reading formulae.......")
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$form_[[1]]<- form[par_names == "lambda0"][[1]]
      private$form_[[2]] <- form[par_names == "sigma"][[1]]
      private$form_[[3]] <- form[par_names == "phi"][[1]]
      private$form_[[4]] <- form[par_names == "beta"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c("lambda0", "sigma", "phi", "beta")
      private$make_par() 
      private$link2response_ <- list(lambda0 = "exp", 
                            sigma = "exp", 
                            phi = "plogis", 
                            beta = "exp", 
                            D = "exp")
      private$response2link_ <- list(lambda0 = "log", 
                                     sigma = "log", 
                                     phi = "qlogis", 
                                     beta = "log", 
                                     D = "log")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$num_cores_ = num_cores
      private$print_ = print 
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
     if (name == "D") {
       D <- do.call(private$link2response_[[5]], list(as.numeric(self$par()["D"])))
       return(D)
     }
     if (name == "beta") {
       k_save <- k
       k <- 2:(private$data_$n_occasions()) 
     } 
     if (name == "phi" & is.null(k)) {
       k <- 1:(private$data_$n_occasions() - 1) 
     }
     covs <- private$data_$covs(j = j, k = k, m = m)
     i_par <- which(names(private$form_) == name) 
     X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
     if (name %in% c("phi", "beta") & "t" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("t1", sep ="")] 
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
      a0 <- self$get_par("beta", j = 1, k = 1)
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1 - a0, a0, 0), nrow = n_mesh, ncol = 3, byrow = TRUE)
      pr0 <- pr0 / n_mesh
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
      tpms <- vector("list", length = n_occasions)
      dt <- diff(private$data_$time())
      for (k in 1:(n_occasions - 1)) {
        phi <- self$get_par("phi", j = 1, k = k, m = 1)
        phi <- phi^dt[k]
        tpms[[k]] <- matrix(c(1 - pr_entry[k], 0, 0, 
                            pr_entry[k], phi, 0,
                            0, 1 - phi, 1), nrow = 3, ncol = 3)
      }
      return(tpms)
    }, 
    
    calc_pr_capture = function() {
      dist <- t(private$data_$distances())
      n_occasions <- private$data_$n_occasions("all")
      n_primary <- private$data_$n_primary()
      S <- private$data_$n_secondary() 
      if (n_primary == 1) {
        n_primary <- n_occasions
        S <- rep(1, n_occasions)
      }
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
                                3, 
                                self$data()$detector_type(),
                                n_primary, 
                                S)
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
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = 3)
        pr_empty[[j]][, 2] <- exp(-t(trap_usage[, j]) %*% enc_rate[j,,])
      }
      # average over all life histories 
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      pdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms, 3); 
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms,
			private$num_cores_, 3, rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_pdet())
      llk <- llk + self$calc_D_llk()
      #plot(self$par()$beta[-1])
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
    
    fit = function(ini_par = NULL) {
      if (!is.null(ini_par)) self$set_par(ini_par)
      par <- self$par()
      w_par <- private$convert_par2vec(par)
      if (private$print_) cat("Fitting model..........\n")
      t0 <- Sys.time()
      mod <- suppressWarnings(optim(w_par, 
                                    private$calc_negllk,
                   names = names(w_par),
                   hessian = TRUE))
      t1 <- Sys.time()
      difft <- t1 - t0
      if (private$print_) cat("Completed model fitting in", difft, attr(difft, "units"), "\n")
      mle <- mod$par
      code <- mod$convergence 
      if (code > 0) warning("model failed to converge with optim code ", code, "\n")
      if (private$print_ & code == 0) cat("Checking convergence.......converged", code, "\n")
      names(mle) <- names(w_par)
      mle <- private$convert_vec2par(mle)
      #mle <- lapply(mle, function(x) {y <- x; names(y) <- NULL; return(y)})
      self$set_par(mle)
      private$mle_ <- mle
      private$llk_ <- -mod$value
      private$V_ <- solve(mod$hessian)
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
    mle_llk = function() {return(private$llk_)},
  
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
    data_ = NULL,
    form_ = NULL, 
    par_ = NULL, 
    link2response_ = NULL, 
    response2link_ = NULL, 
    mle_ = NULL,
    results_ = NULL,
    D_tab_ = NULL, 
    Dk_ = NULL, 
    var_ = NULL, 
    confint_ = NULL, 
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05,
    print_ = NULL, 
    num_cores_ = NULL, 
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_par <- numeric(4)
      private$par_ <- vector(mode = "list", length = 5)
      for (par in 1:4) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        if (par %in% c(3, 4) & "t" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("t", private$data_$n_occasions() - 1, sep ="")]
        } else {
          n_par[par] <- ncol(X)
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)
        }
        private$par_[[par]] <- par_vec
      }
      private$par_[[5]] <- 0
      names(private$par_) <- c(names(private$form_), "D") 
    }, 
    
    initialise_par = function(start) {
        private$par_$lambda0[1] <- do.call(private$response2link_$lambda0, 
                                           list(start$lambda0))
        private$par_$sigma[1] <-do.call(private$response2link_$sigma, 
                                           list(start$sigma))
        private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
        private$par_$beta[1] <- do.call(private$response2link_$beta,
                                        list(c(1 / start$beta - 1)/self$data()$n_occasions()))
        private$par_$D <- do.call(private$response2link_$D, 
                                           list(start$D))
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
      print(results)
      colnames(results) <- c("Estimate", "SE", "LCL", "UCL")
      rownames(results) <- names(par)
      private$results_ <- results 
      
      ## D tab
      D_tab <- matrix(0, nr = self$data()$n_occasions(), nc = 4)
      colnames(D_tab) <- c("Estimate", "SE", "LCL", "UCL")
      rownames(D_tab) <- private$data_$time()
      D_tab[, 1] <- private$Dk_
      D_tab[ ,2] <- sqrt((exp(private$var_$Dkvar) - 1) * private$Dk_^2 * exp(private$var_$Dkvar))
      D_tab[, 3] <- private$confint_$Dk$lcl
      D_tab[, 4] <- private$confint_$Dk$ucl

      private$D_tab_ <- D_tab
    
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
    del_pdet <- grad(private$calc_wpdet, wpar[-np]) 
    # get covariance matrix 
    V <- self$cov_matrix()
    sds <- sqrt(diag(V))
    # theta variance
    dist <- self$data()$distances()
    V_theta <- V[-np, -np] * pdet 
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
        del_alpha <- grad(private$calc_alpha, wpar[-np], k = k)
        sig_alpha <- t(del_alpha) %*% V_theta %*% del_alpha
        Dk_var[k] <- Dk[k] * Dvar / self$get_par("D") + sig_alpha / self$get_par("D")
      }
      Dk_var <- Dk_var / (self$data()$area() * Dk)
    }
    private$var_ <- list(sds = sds, Dvar = Dvar, Dkvar = Dk_var)
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
        varlog <- 1 + log(private$var_$Dkvar / Dk^2)
        Dk_C <- exp(qnorm(lev) * sqrt(private$var_$Dkvar))
        Dk_lcl <- Dk / Dk_C
        Dk_ucl <- Dk * Dk_C
        Dkmat <- data.frame(Dk = Dk, lcl = Dk_lcl, ucl = Dk_ucl) 
      }
      private$confint_ <- list(est = estmat, Dk = Dkmat)
    },
  
     infer_D = function() {
       tpms <- self$calc_tpms()
       a0 <- self$get_par("beta", k = 1)
       pr0 <- c(1 - a0, a0, 0)
       private$Dk_ <- C_calc_D(self$get_par("D"), self$data()$n_occasions(), pr0, tpms)
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
      par$phi <- vec[grep("phi", names)]
      names(par$phi) <- gsub("phi.", "", names(par$phi))
      par$beta <- vec[grep("beta", names)]
      names(par$beta) <- gsub("beta.", "", names(par$beta))
      par$D <- vec["D"]
      names(par$D) <- NULL 
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



