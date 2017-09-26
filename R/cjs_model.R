# Copyright (c) 2017 Richard Glennie
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
# openpopscr project: open population spatial capture-recapture 
#
# file description: Cormack-Jolly-Seber model class
#
################################################################################

#' Cormack-Jolly-Seber model class 
#' 
#' @description Encapsulates Cormack-Jolly-Seber model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item scr_data: a ScrData object 
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item start: a named list of starting values 
#'   \item num_cores (optional, default = 1): number of processors cores to use 
#'   in parallelised code 
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k at mesh point m (if j, k, or m omitted returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses. Note, the model will simulate 
#'    data using this parameter, but will only present inference based on the maximum likelihood
#'    estimates. 
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_pr_entry(): computes vector with entry j equal to probability of individual unborn up to occasion j
#'  being born just after occasion j
#'  \item calc_tpms(): returns list of transition probability matrix for each occasion 
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit(nsims = 99): fit the model by obtaining the maximum likelihood estimates. Estimates of
#'        density are obtained from parametric boostrap with nsim resamples. 
#'  \item simulate(): simulate ScrData object from fitted model
#'  \item par(): return current parameter of the model 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#'  \item sample_D(nsims): generate nsims time series for density, sampled from fitted model
#'  \item sample_R(nsims): generate nsims time series of arrival numbers, sampled from fitted model  
#' }
#' 
CjsModel <- R6Class("JsModel", 
  public = list(
    
    initialize = function(form, data, start, num_cores = 1) {
      private$data_ <- data
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$form_[[1]]<- form[par_names == "lambda0"][[1]]
      private$form_[[2]] <- form[par_names == "sigma"][[1]]
      private$form_[[3]] <- form[par_names == "phi"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c("lambda0", "sigma", "phi")
      private$make_par() 
      private$link2response_ <- list(lambda0 = "exp", 
                            sigma = "exp", 
                            phi = "plogis")
      private$response2link_ <- list(lambda0 = "log", 
                                     sigma = "log", 
                                     phi = "qlogis") 
      private$initialise_par(start)
      private$num_cores_ = num_cores
    },
    
    get_par = function(name, j = NULL, k = NULL, m = NULL) {
     if (name == "phi" & is.null(k)) {
       k <- 1:(private$data_$n_occasions() - 1) 
     }
     covs <- private$data_$covs(j = j, k = k, m = m)
     i_par <- which(names(private$form_) == name) 
     X <- model.matrix(private$form_[[i_par]], data = as.data.frame(covs)) 
     if (name %in% c("phi") & "t" %in% all.vars(private$form_[[i_par]])) {
       X <- X[,colnames(X) != paste("t1", sep ="")] 
     }
     theta <- private$par_[[i_par]]
     l_par <- which(names(private$link2response_) == name)
     resp <- do.call(private$link2response_[[l_par]], list(X %*% theta))
     return(resp)
    }, 
    
    set_par = function(par) {
      private$par_ <- par
    },
    
    set_mle = function(mle, V, llk, nsims = 999) {
      private$mle_ <- mle
      private$llk_ <- -llk
      private$V_ <- V
      sd <- sqrt(diag(private$V_))
      names(sd) <- names(private$convert_par2vec(private$convert_natural2working(mle)))
      sd <- private$convert_working2natural(private$convert_vec2par(sd))
      confints <- private$calc_confint()
      private$make_results(sd, confints)
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
      tpms <- vector("list", length = n_occasions)
      dt <- diff(private$data_$time())
      for (k in 1:(n_occasions - 1)) {
        phi <- self$get_par("phi", j = 1, k = k, m = 1)
        if (dt[k] > 1) phi <- phi^dt[k]
        tpms[[k]] <- matrix(c(phi, 0, 
                            0, 1 - phi), nrow = 2, ncol = 2)
      }
      return(tpms)
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
      prob <- C_calc_pr_capture(n, n_occasions, n_traps, n_meshpts, capthist, 
                               enc_rate0, trap_usage, private$num_cores_, 2)
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
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = 2)
        pr_empty[[j]][, 1] <- exp(-t(trap_usage[, j]) %*% enc_rate[j,,])
      }
      # average over all life histories 
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      pdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms); 
      return(pdet)
    },
    
    calc_llk = function(param = NULL, names = NULL) {
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, private$num_cores_, 2)
      cat("llk:", llk, "\n")
      return(llk)
    },
    
    fit = function(ini_par = NULL) {
      if (!is.null(ini_par)) self$set_par(ini_par)
      par <- self$par()
      w_par <- private$convert_par2vec(par)
      mod <- optim(w_par,
                   private$calc_negllk,
                   names = names(w_par),
                   hessian = TRUE)
      mle <- mod$par
      names(mle) <- names(w_par)
      mle <- private$convert_vec2par(mle)
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
    }
    options(scipen = 0)
  }, 
  
  simulate = function(n = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (!is.null(n)) n <- private$data_$n()
    num.meshpts <- private$data_$n_meshpts()
    mesh <- private$data_$mesh()
    D <- n * private$data_$area() * 100 
    # simulate population
    pop <- sim.popn(D = D, core = mesh, Ndist = "fixed", buffertype = "convex")
    n_occasions <- private$data_$n_occasions()
    dt <- diff(private$data_$time())
    phi <- self$get_par("phi", j = 1, m = 1)
    phi <- phi ^ dt
    life <- matrix(0, nr = nrow(pop), ncol = n_occasions) 
    surv <- rgeom(n, 1 - phi)
    for (i in 1:n) life[i, 1:(1 + surv)] <- 1   
    lambda0 <- self$get_par("lambda0", m = 1)
    sigma <- self$get_par("sigma", j = 1, m = 1)
    capture_history <- sim.capthist(private$data_$traps(), 
                                    popn = pop, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = lambda0, 
                                                     sigma = sigma), 
                                    noccasions = n_occasions,
                                    renumber = FALSE)
    # thin capture history by survival
    n <- dim(capture_history)[1]
    ids <- as.numeric(rownames(capture_history))
    life <- life[ids,]
    birth_time <- birth_time[ids]
    seen <- rep(TRUE, n)
    for (i in seq(n)) {
      capture_history[i, ,] <- diag(life[i,]) %*% capture_history[i, ,]
    } 
    new_dat <- ScrData$new(capture_history, mesh, private$data_$time())
    return(new_dat)
  }, 
  
  ## UP TO HERE 
  
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
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05, 
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
                                        list(c(1 / start$beta - 1)))
        if (length(private$par_$beta) > 1) {
          private$par_$beta <- rep(1 / start$beta - 1, private$data_$n_occasions() - 1)
          private$par_$beta <- private$par_$beta / length(private$par_$beta)
          private$par_$beta <- do.call(private$response2link_$beta, 
                                           list(private$par_$beta))
          private$par_$beta[-1] <- private$par_$beta[-1] - private$par_$beta[1]
        }
        private$par_$D <- do.call(private$response2link_$D, 
                                           list(start$D))
    }, 
    
    make_results = function(nsims) {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      par <- private$convert_par2vec(private$mle_) 
      sd <- sqrt(diag(private$V_))
      ci <- private$calc_confint()
      Dinfer <- private$infer_D(nsims = nsims)
      results <- cbind(par, sd, ci$LCL, ci$UCL)
      colnames(results) <- c("Estimate", "Std. Error", "LCL", "UCL")
      rownames(results) <- names(par)
      private$results_ <- results 
      
      ## D tab
      D_tab <- matrix(0, nr = private$data_$n_occasions(), nc = 3)
      colnames(D_tab) <- c("Estimate", "LCL", "UCL")
      rownames(D_tab) <- private$data_$time()
      D_tab[, 1] <- Dinfer$mean
      D_tab[, 2] <- Dinfer$lcl
      D_tab[, 3] <- Dinfer$ucl

      private$D_tab_ <- D_tab
    
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
  
     infer_D = function(nsims, extract_samples = 0) {
       save_par <- self$par()
       self$set_par(private$mle_)
       sd <- sqrt(diag(private$V_))
       w_par <- private$convert_par2vec(self$par())
       names <- names(w_par)
       tpms <- self$calc_tpms()
       a0 <- self$get_par("beta", j = 1, k = 1, m = 1)
       pr0 <- c(1 - a0, a0, 0)
       n_occasions <- private$data_$n_occasions()
       D <- do.call(private$link2response_$D, list(self$par()$D))
       mean <- C_calc_D(D, n_occasions, pr0, tpms)
       Dsims <- matrix(0, nr = nsims, nc = n_occasions)
       Rsims <- matrix(0, nr = nsims, nc = n_occasions)
       for (sim in 1:nsims) {
         new_par <- rnorm(length(w_par), w_par, sqrt(diag(private$V_)))
         names(new_par) <- names
         self$set_par(private$convert_vec2par(new_par))
         tpms <- self$calc_tpms()
         a0 <- self$get_par("beta", j = 1, k = 1, m = 1)
         pr0 <- c(1 - a0, a0, 0)
         D <- do.call(private$link2response_$D, list(self$par()$D))
         Dsims[sim, ] <- C_calc_D(D, n_occasions, pr0, tpms)
         Rsims[sim, ] <- self$get_par("beta", m = 1) * D
       }
       if (extract_samples == 1) return(Dsims)
       if (extract_samples ==  2) return(Rsims)  
       var <- apply(Dsims, 2, var)
       alp <- qnorm(1 - private$sig_level_ / 2)
       sd <- sqrt(var)
       lcl <- apply(Dsims, 2, function(x){return(quantile(x, private$sig_level_ / 2))})
       ucl <- apply(Dsims, 2, function(x){return(quantile(x, 1 - private$sig_level_ / 2))})
       
       self$set_par(save_par)
       return(list(mean = mean, 
                   sd = sd, 
                   lcl = lcl, 
                   ucl = ucl))
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



