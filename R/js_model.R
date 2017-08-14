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
# file description: Jolly-Seber model class
#
################################################################################

#' Jolly-Seber model class 
#' 
#' @description Encapsulates Jolly-Seber model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item scr_data: a ScrData object 
#'   \item par: a named list of initial parameter values
#'   \item num_cores (optional, default = 1): number of processors cores to use 
#'   in parallelised code 
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item set_par(): can change the parameter the model uses. Note, the model will simulate 
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
#'  \item fit(): fit the model by obtaining the maximum likelihood estimates
#'  \item plot(): plot estimated density over time with 95% confidence intervals
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
    
    initialize = function(scr_data, par, num_cores = 1) {
      private$data_ = scr_data
      private$par_ = par
      private$mle_ = NULL 
      private$num_cores_ = num_cores
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
      Dinference <- private$infer_D(nsims)
      private$make_results(sd, confints, Dinference)
    }, 
    
    calc_D_llk = function() {
      D <- self$par()$D
      A <- private$data_$area()
      n <- private$data_$n()
      pdet <- self$calc_pdet()
      llk <- n * log(D * A * pdet) - D * A * pdet - lfactorial(n)
      names(llk) <- NULL 
      return(llk)
    },
    
    calc_initial_distribution = function() {
      a0 <- self$par()$beta[1]
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1 - a0, a0, 0), nrow = n_mesh, ncol = 3, byrow = TRUE)
      pr0 <- pr0 / n_mesh
      return(pr0)
    },
    
    calc_pr_entry = function() {
      beta <- self$par()$beta 
      n_occasions <- private$data_$n_occasions()  
      pr_entry <- rep(0, n_occasions - 1)
      prod <- 1 - beta[1]
      if (length(beta) != n_occasions) {
        dt <- diff(private$data_$time())
        beta <- c(beta[1], (1 - beta[1]) * dt / sum(dt))
      }
      for (j in 1:(n_occasions - 1)) {
        pr_entry[j] <- beta[j + 1] / prod
        prod <- prod * (1 - pr_entry[j])
      }
      return(pr_entry)
    }, 
    
    calc_tpms = function() {
      # compute entry probabilities 
      pr_entry <- self$calc_pr_entry()      
      n_occasions <- private$data_$n_occasions()
      tpms <- vector("list", length = n_occasions)
      phi <- self$par()$phi
      if (length(phi) == 1) {
        phi <- rep(phi, private$data_$n_occasions() - 1)
        if (!is.null(private$data_$time())) phi <- phi^diff(private$data_$time())
      } 
      for (j in 1:(n_occasions - 1)) {
        tpms[[j]] <- matrix(c(1 - pr_entry[j], 0, 0, 
                            pr_entry[j], phi[j], 0,
                            0, 1 - phi[j], 1), nrow = 3, ncol = 3)
      }
      return(tpms)
    }, 
    
    calc_pr_capture = function() {
      dist <- t(private$data_$distances())
      n_occasions <- private$data_$n_occasions()
      lambda0 <- self$par()$lambda0
      if (length(lambda0) == 1) lambda0 <- rep(lambda0, n_occasions)
      sigma <- self$par()$sigma
      if (length(sigma) == 1) sigma <- rep(sigma, n_occasions)
      enc_rate0 <- array(0, dim = c(nrow(dist), ncol(dist), n_occasions)) 
      for (j in 1:n_occasions) {
        enc_rate0[,,j] <- lambda0[j] * exp(-dist ^ 2 / (2  * sigma[j] ^ 2))
      }
      trap_usage <- usage(private$data_$traps())
      n <- private$data_$n()
      n_meshpts <- private$data_$n_meshpts() 
      n_traps <- private$data_$n_traps()
      capthist <- private$data_$capthist()
      prob <- C_calc_pr_capture(n, n_occasions, n_traps, n_meshpts, capthist, 
                               enc_rate0, trap_usage, private$num_cores_)
      return(prob)
    },
    
    calc_pdet = function() {
      # compute probability of zero capture history 
      dist <- private$data_$distances()
      n_occasions <- private$data_$n_occasions()
      lambda0 <- self$par()$lambda0
      if (length(lambda0) == 1) lambda0 <- rep(lambda0, n_occasions)
      sigma <- self$par()$sigma
      if (length(sigma) == 1) sigma <- rep(sigma, n_occasions)
      enc_rate <- array(0, dim = c(n_occasions, nrow(dist), ncol(dist))) 
      for (j in 1:n_occasions) {
        enc_rate[j,,] <- lambda0[j] * exp(-dist ^ 2 / (2  * sigma[j] ^ 2))
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
      pdet <- C_calc_pdet(private$data_$n_occasions(), pr0, pr_empty, tpms); 
      return(pdet)
    },
    
    calc_llk = function(param = NULL, names = NULL, scale = "natural") {
      if (!is.null(param)) {
        if (scale == "natural") {
          if(!is.null(names)) names(param) <- names 
          self$set_par(param)  
        } else {
          if (scale == "working") {
            if(!is.null(names)) names(param) <- names 
            par <- private$convert_vec2par(param)
            self$set_par(private$convert_working2natural(par))
          }
        }
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
      llk <- C_calc_llk(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, private$num_cores_)
      # compute log-likelihood
      llk <- llk - n * log(self$calc_pdet())
      llk <- llk + self$calc_D_llk()
      #plot(self$par()$beta[-1])
      cat("llk:", llk, 
          "lambda0:", mean(self$par()$lambda0), 
          "sigma:", mean(self$par()$sigma), 
          "phi:", mean(self$par()$phi),
          "beta:", self$par()$beta[1], 
          "D :", self$par()$D, "\n")
      return(llk)
    },
    
    fit = function(ini_par = NULL, nsims = 999) {
      if (!is.null(ini_par)) self$set_par(ini_par)
      par <- self$par()
      w_par <- private$convert_natural2working(par)
      w_par <- private$convert_par2vec(w_par)
      # mod <- nlm(private$calc_negllk,
      #            w_par,
      #            names = names(w_par),
      #            scale = "working",
      #            hessian = TRUE)
      mod <- optim(w_par,
                   private$calc_negllk,
                   names = names(w_par),
                   scale = "working",
                   hessian = TRUE)
      mle <- mod$par
      #mle <- mod$estimate
      #mle <- w_par
      names(mle) <- names(w_par)
      mle <- private$convert_vec2par(mle)
      mle <- private$convert_working2natural(mle)
      mle <- lapply(mle, function(x) {y <- x; names(y) <- NULL; return(y)})
      private$mle_ <- mle
      private$llk_ <- -mod$value
      private$V_ <- solve(mod$hessian)
      if (any(diag(private$V) <= 0)) {
        cat("Variance estimates not reliable, do a bootstrap.")
        return(0)
      } else {
      #private$V_ <- matrix(runif(length(w_par)^2), nr = length(w_par))
        sd <- sqrt(diag(private$V_))
        names(sd) <- names(w_par)
        sd <- private$convert_working2natural(private$convert_vec2par(sd))
      #private$llk_ <- -mod$minimum
        confints <- private$calc_confint()
        Dinference <- private$infer_D(nsims)
        private$make_results(sd, confints, Dinference)
      }
    }, 
    
  print = function() {
    options(scipen = 999)
    if (is.null(private$mle_)) {
      print("Fit model using $fit method")
    } else {
      cat("DETECTION \n")
      print(signif(private$tab_, 4))
    
      cat("--------------------------------------------------------------------------------")
      cat("\n SURVIVAL \n")
      print(signif(private$phi_tab_, 4))
      
      cat("--------------------------------------------------------------------------------")
      cat("\n ARRIVAL \n")
      if (!is.null(private$beta_tab_)) {
        print(signif(private$beta_tab_, 4))
      } else {
        cat("\n Assumed constant. \n")
      }
      cat("--------------------------------------------------------------------------------")
      cat("\n DENSITY \n")
      print(signif(private$D_tab_, 4))
      cat("--------------------------------------------------------------------------------")
    }
    options(scipen = 0)
  }, 
  
  plot = function(par = "D") {
    
    plot_var <- function(dat) {
      time <- 1:private$data_$n_occasions()
       if (!is.null(private$data_$time())) time <- private$data_$time()
       full_time <- seq(time[1], time[length(time)], 1)
       dat <- data.frame(mean = dat[, 1], 
                           lcl = dat[, 2], 
                           ucl = dat[ ,3], 
                           t = time)
       grp <- ggplot2::ggplot(dat, ggplot2::aes(x = t)) + 
         ggplot2::geom_line(ggplot2::aes(y = mean), col = "black") + 
         ggplot2::geom_point(ggplot2::aes(y = mean), col = "black") +
         ggplot2::geom_ribbon(ggplot2::aes(ymin = lcl, ymax = ucl), linetype = 2, alpha = 0.2) + 
         ggplot2::scale_x_continuous("Occasion", breaks = full_time, labels = full_time) + 
         ggplot2::scale_y_continuous("Density (per square kilometre)") + 
         ggplot2::theme_bw() +
         ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                        panel.grid.minor = ggplot2::element_blank())
       grp
    }
    
   if (is.null(private$mle_)) {
     private$data_$plot()
   } else {
     if (par == "D") {
       plot_var(private$D_tab_)
       
     } else {
       if (is.null(self$par()[par])) stop("parameter does not exist")
       if (length(self$par()[par]) == 1) stop("parameter not time-varying")
       plot_var(self$par()[par])
     }
   }
  }, 
 
  simulate = function(seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    num.meshpts <- private$data_$n_meshpts()
    mesh <- private$data_$mesh()
    D <- self$par()$D / 100
    # simulate population
    pop <- sim.popn(D = D, core = mesh, Ndist = "poisson", buffertype = "convex")
    a0 <- self$par()$beta[1]
    birth_time <- sample(c(0, 1), size = nrow(pop), prob = c(1 - a0, a0), replace = TRUE) 
    beta <- self$par()$beta[-1]
    n_occasions <- private$data_$n_occasions()
    if (length(beta) == 0) {
      beta <- rep(1 , n_occasions - 1) 
      if (!is.null(private$data_$time())) beta <- diff(private$data_$time())
      beta <- beta / sum(beta) 
    }
    birth_time <- ifelse(birth_time == 0, 
                         sample(2:n_occasions, 
                                size = nrow(pop), 
                                prob = beta,
                                replace = TRUE), 
                         1)
    
    phi <- self$par()$phi 
    if (length(phi) == 1) {
      phi <- rep(phi, n_occasions - 1)
      if (!is.null(private$data_$time())) phi <- phi^diff(private$data_$time())
    }
    life <- matrix(0, nr = nrow(pop), ncol = n_occasions) 
    life[birth_time == 1, 1] <- 1
    for (j in 2:n_occasions) {
      life[birth_time < j,j] <- rbinom(sum(birth_time < j), 1, phi[j - 1]) 
      life[birth_time == j, j] <- 1
    }
    # generate capture histories
    plot(mesh)
    plot(pop[life[,1] == 1, ], add = T)
    plot(private$data_$traps(), add = T)
    capture_history <- sim.capthist(private$data_$traps(), 
                                    popn = pop, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = self$par()$lambda0, 
                                                     sigma = self$par()$sigma), 
                                    noccasions = n_occasions,
                                    renumber = FALSE)
    # thin capture history by survival
    n <- dim(capture_history)[1]
    ids <- as.numeric(rownames(capture_history))
    life <- life[ids,]
    birth_time <- birth_time[ids]
    seen <- rep(TRUE, n)
    for (i in seq(n)) {
      life[i, birth_time[i]:n_occasions] <- cumprod(life[i, birth_time[i]:n_occasions])
      capture_history[i, ,] <- diag(life[i,]) %*% capture_history[i, ,]
      if (sum(capture_history[i, ,]) == 0) seen[i] <- FALSE
    } 
    capture_history <- subset(capture_history, (1:n)[seen])
    new_dat <- ScrData$new(capture_history, mesh, private$data_$time())
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
        ests$estimates <- private$tab_ 
        ests$phi <- private$phi_tab_ 
        if (!is.null(private$beta_tab_)) ests$beta <- private$beta_tab_
        ests$D <- private$D_tab_
      }
      return(ests)
    },
    cov_matrix = function() {return(private$V_)}, 
    mle_llk = function() {return(private$llk_)},
    sample_D = function(nsims = 999) {
      if (is.null(private$mle_)) stop("Fit model using $fit method.")
      sds <- sqrt(diag(private$V_))
      return(private$infer_D(nsims, sds, extract_samples = 1))
    }, 
     
    sample_R = function(nsims = 999) {
      if (is.null(private$mle_)) stop("Fit model using $fit method.")
      sds <- sqrt(diag(private$V_))
      return(private$infer_D(nsims, sds, extract_samples = 2))
    } 
  
),
                   
  private = list(
    data_ = NULL,
    par_ = NULL, 
    mle_ = NULL,
    tab_ = NULL,
    phi_tab_ = NULL, 
    beta_tab_ = NULL, 
    D_tab_ = NULL, 
    V_ = NULL, 
    llk_ = NULL, 
    sig_level_ = 0.05, 
    num_cores_ = NULL, 
    
    make_results = function(sd, confint, Dinference) {
      if (is.null(private$mle_)) print("Fit model using $fit method.")
      n_occasions <- private$data_$n_occasions()
      tab_n_rows <- n_occasions * 2
      shift <- rep(n_occasions, 2)
      if (length(private$mle_$lambda0) == 1) {
        tab_n_rows <- tab_n_rows - n_occasions + 1
        shift[1] <- 1 
      }
      if (length(private$mle_$sigma) == 1) {
        tab_n_rows <- tab_n_rows - n_occasions + 1
        shift[2] <- 1 
      }
      tab <- matrix(0, nr = tab_n_rows, nc = 3)
      colnames(tab) <- c("Est.", "LCL", "UCL")
      tab[1:shift[1], ] <- c(private$mle_$lambda0, 
                             confint$LCL$lambda0, 
                             confint$UCL$lambda0)
      
      tab[(shift[1] + 1):(shift[1] + shift[2]), ] <- c(private$mle_$sigma, 
                                                       confint$LCL$sigma, 
                                                       confint$UCL$sigma)
      
      row_names <- NULL 
      if (length(private$mle_$lambda0) == 1) {
        row_names <- c(row_names, "lambda0")
      } else {
        row_names <- c(paste("lambda0_", 1:n_occasions, sep = ""))
      }
      
      if (length(private$mle_$sigma) == 1) {
        row_names <- c(row_names, "sigma")
      } else {
        row_names <- c(row_names, paste("sigma_", 1:n_occasions, sep = ""))
      }
      
      rownames(tab) <- row_names 
      
      
      if (length(self$par()$phi) == 1) {
        phi_tab <- matrix(0, nr = 1, nc = 3)
      } else {
        phi_tab <- matrix(0, nr = private$data_$n_occasions() - 1, nc = 3)
        rownames(phi_tab) <- 1:(private$data_$n_occasions() - 1)
      }
      colnames(phi_tab) <- c("Est.", "LCL", "UCL")
      phi_tab[, 1] <- private$mle_$phi
      phi_tab[, 2] <- confint$LCL$phi
      phi_tab[, 3] <- confint$UCL$phi 
      
      beta_tab <- NULL 
      beta_tab <- matrix(0, nr = length(self$par()$beta), nc = 3)
      colnames(beta_tab) <- c("Est.", "LCL", "UCL")
      rownames(beta_tab) <- 1:(length(self$par()$beta))
      beta_tab[, 1] <- private$mle_$beta
      beta_tab[, 2] <- Dinference$beta_lcl
      beta_tab[, 3] <- Dinference$beta_ucl 
      
      D_tab <- matrix(0, nr = private$data_$n_occasions(), nc = 3)
      colnames(D_tab) <- c("Est.", "LCL", "UCL")
      rownames(D_tab) <- 1:(private$data_$n_occasions())
      if (!is.null(private$data_$time())) rownames(D_tab) <- private$data_$time()
      D_tab[, 1] <- Dinference$mean
      D_tab[, 2] <- Dinference$lcl
      D_tab[, 3] <- Dinference$ucl
      
      # save results 
      private$tab_ <- tab
      private$phi_tab_ <- phi_tab
      private$beta_tab_ <- beta_tab
      private$D_tab_ <- D_tab
    
  }, 
  
   calc_confint = function() {
      V <- private$V_ 
      sds <- sqrt(diag(V))
      est <- private$convert_natural2working(private$mle_)
      est <- private$convert_par2vec(est)
      alp <- qnorm(1 - private$sig_level_ / 2)
      lcl <- est - alp * sds
      ucl <- est + alp * sds 
      names(lcl) <- names(ucl) <- names(est)
      lcl <- private$convert_working2natural(private$convert_vec2par(lcl))
      ucl <- private$convert_working2natural(private$convert_vec2par(ucl))
      return(list(LCL = lcl, UCL = ucl))
    },
  
     calc_Dt = function(wpar, J) {
       par <- private$convert_working2natural(private$convert_vec2par(wpar))
       tpms <- self$calc_tpms()
       pr0 <- c(1 - par$beta[1], par$beta[1], 0)
       return(C_calc_Dt(self$par()$D, J, pr0, tpms))
     },
  
     infer_D = function(nsims, sd, extract_samples = 0) {
       save_par <- self$par()
       self$set_par(private$mle_)
       w_par <- private$convert_par2vec(private$convert_natural2working(self$par()))
       names <- names(w_par)
       # compute transition probability matrices 
       tpms <- self$calc_tpms()
       # initial distribution 
       #pr0 <- self$calc_initial_distribution()
       pr0 <- c(1 - self$par()$beta[1], self$par()$beta[1], 0)
       n_occasions <- private$data_$n_occasions()
       mean <- C_calc_D(self$par()$D, n_occasions, pr0, tpms)
       Dsims <- matrix(0, nr = nsims, nc = n_occasions)
       Rsims <- matrix(0, nr = nsims, nc = n_occasions)
       beta_sims <- matrix(0, nr = nsims, nc = length(self$par()$beta))
       for (sim in 1:nsims) {
         new_par <- rnorm(length(w_par), w_par, sqrt(diag(private$V_)))
         names(new_par) <- names
         self$set_par(private$convert_working2natural(private$convert_vec2par(new_par)))
         beta_sims[sim, ] <- self$par()$beta
         tpms <- self$calc_tpms()
         pr0 <- c(1 - self$par()$beta[1], self$par()$beta[1], 0)
         Dsims[sim, ] <- C_calc_D(self$par()$D, n_occasions, pr0, tpms)
         Rsims[sim, ] <- beta_sims[sim, ] * self$par()$D
       }
       if (extract_samples == 1) return(Dsims)
       if (extract_samples ==  2) return(Rsims)  
       var <- apply(Dsims, 2, var)
       alp <- qnorm(1 - private$sig_level_ / 2)
       sd <- sqrt(var)
       lcl <- apply(Dsims, 2, function(x){return(quantile(x, private$sig_level_ / 2))})
       ucl <- apply(Dsims, 2, function(x){return(quantile(x, 1 - private$sig_level_ / 2))})
       
       beta_var <- apply(beta_sims, 2, var)
       beta_sd <- sqrt(beta_var)
       beta_lcl <- apply(beta_sims, 2, function(x){return(quantile(x, private$sig_level_ / 2))})
       beta_ucl <- apply(beta_sims, 2, function(x){return(quantile(x, 1 - private$sig_level_ / 2))})
       
       self$set_par(save_par)
       return(list(mean = mean, 
                   sd = sd, 
                   lcl = lcl, 
                   ucl = ucl,
                   beta_sd = beta_sd, 
                   beta_lcl = beta_lcl,
                   beta_ucl = beta_ucl))
    },
  
     calc_negllk = function(param = NULL, names = NULL, scale = "natural") {
       # if(!is.null(names)) names(param) <- names 
       # par <- private$convert_vec2par(param)
       # self$set_par(private$convert_working2natural(par))
       negllk <- -self$calc_llk(param, names, scale)
       # if (!is.null(self$par()$beta)) {
       #   n_param <- length(param)
       #   lagrange <- param[n_param]
       #   negllk <- negllk + lagrange * (1 - sum(self$par()$beta))  
       # }
       return(negllk)
    },
  
    convert_working2natural = function(wpar) {
       par <- NULL
       if (length(wpar$beta) == 1) {
         par$beta <- plogis(wpar$beta)
       } else {
         exp_beta <- exp(wpar$beta)
         par$beta <- c(1, exp_beta) / (1 + sum(exp_beta))
       }
       par$sigma <- exp(wpar$sigma)
       par$lambda0 <- exp(wpar$lambda0)
       par$phi <- plogis(wpar$phi)
       par$D <- exp(wpar$D)
       return(par)
    }, 
    
    convert_natural2working = function(par) {
      wpar <- NULL
      if (length(par$beta) == 1) {
        wpar$beta <- private$qlog(par$beta[1])
      } else {
        wpar$beta <- log(par$beta[-1] / par$beta[1])
      }
      wpar$sigma <- log(par$sigma)
      wpar$lambda0 <- log(par$lambda0)
      wpar$phi <- private$qlog(par$phi)
      wpar$D <- log(par$D)
      return(wpar)
    },
  
   convert_par2vec = function(par) {
      return(unlist(par))
    },
    
    convert_vec2par = function(vec) {
      par <- NULL
      n_occasions <- private$data_$n_occasions()
      names <- names(vec)
      if ("beta1" %in% names) par$beta <- vec[paste("beta", seq(1, n_occasions - 1), sep = "")]
      if ("beta" %in% names) par$beta <- vec["beta"]
      if ("phi1" %in% names) par$phi <- vec[paste("phi", seq(1, n_occasions - 1), sep ="")]
      if ("lambda01" %in% names) par$lambda0 <- vec[paste("lambda0", seq(1, n_occasions), sep ="")]
      if ("sigma1" %in% names) par$sigma <- vec[paste("sigma", seq(1, n_occasions), sep ="")]
      if ("phi" %in% names) par$phi <- vec["phi"]
      if ("sigma" %in% names) par$sigma <- vec["sigma"]
      if ("lambda0" %in% names) par$lambda0 <- vec["lambda0"]
      par$D <- vec["D"]
      par <- lapply(par, function(v) {w <- v; names(w) <- NULL; w})
      return(par)
    }, 
  
   qlog = function(p) {
     q <- p
     q <- ifelse(abs(1 - q) < 1e-16, 1 - 1e-16, q)
     q <- ifelse(abs(q) < 1e-16, 1e-16, q)
     return(qlogis(q))
   },
    
   on_par_boundary = function() {
     tol <- 1e-10
     mle <- private$mle_
     side <- NULL 
     if (length(mle$beta) > 1) {
       side$beta <- ifelse(mle$beta[-1] < tol, 1, ifelse(mle$beta[-1] > 1 - tol, -1, NA))
     } else {
       side$beta <- NA 
       if (mle$beta < tol) side$beta <- 1
       if (mle$beta > 1 - tol) side$beta <- -1 
     }
     side$phi <- ifelse(mle$phi < tol, 1, ifelse(mle$phi > 1 - tol, -1, NA))
     side$sigma <- NA
     side$sigma <- ifelse(mle$sigma < tol, 1, NA)
     side$lambda0 <- ifelse(mle$lambda < tol, 1, NA)
     side$D <- ifelse(mle$D < tol, 1, NA)
     return(private$convert_par2vec(side))
   }
  
  )                 
)



