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

#' Jolly-Seber with transience model class 
#' 
#' @description Jolly-Seber model with transience: fits model, formats inference, and 
#' simulates from fitted model. This model inherits all the functions of the 
#' JsModel class. Functions that are different or additional are documented here. 
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
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit(nsims = 99): fit the model by obtaining the maximum likelihood estimates. Estimates of
#'        density are obtained from parametric boostrap with nsim resamples. 
#'  \item simulate(): simulate ScrData object from fitted model
#' }
#' 
JsTransientModel <- R6Class("JsTransientModel", 
                            inherit = JsModel, 
  public = list(
    
    initialize = function(form, data, start, num_cores = 1, print = TRUE) {
      if (print) cat("Creating rectangular mesh......")
      newmesh <- rectangularMask(data$mesh())
      private$dx_ <- attr(newmesh, "spacing")
      private$inside_ <- as.numeric(pointsInPolygon(newmesh, data$mesh()))
      cov_list <- data$get_cov_list() 
      private$data_ <- ScrData$new(data$capthist(), newmesh, data$time(), cov_list$cov, cov_list$cov_type)
      box <- attributes(newmesh)$boundingbox
      region <- c(diff(box[1:2, 1]), diff(box[c(1, 3), 2]))
      private$num_cells_ <- numeric(3)
      private$num_cells_[1] <- nrow(newmesh)
      private$num_cells_[2] <- floor(region[1] / private$dx_)
      private$num_cells_[3] <- nrow(newmesh) / private$num_cells_[2]
      if (print) cat("done\n")
      if (print) cat("Reading formulae.......")
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$form_[[1]]<- form[par_names == "lambda0"][[1]]
      private$form_[[2]] <- form[par_names == "sigma"][[1]]
      private$form_[[3]] <- form[par_names == "phi"][[1]]
      private$form_[[4]] <- form[par_names == "beta"][[1]]
      private$form_[[5]] <- form[par_names == "sd"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c("lambda0", "sigma", "phi", "beta", "sd")
      private$make_par() 
      private$link2response_ <- list(lambda0 = "exp", 
                            sigma = "exp", 
                            phi = "plogis", 
                            beta = "exp", 
                            sd = "exp", 
                            D = "exp")
      private$response2link_ <- list(lambda0 = "log", 
                                     sigma = "log", 
                                     phi = "qlogis", 
                                     beta = "log",
                                     sd = "log", 
                                     D = "log")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$num_cores_ = num_cores
      private$print_ = print 
    },
    
    calc_initial_distribution = function() {
      a0 <- self$get_par("beta", j = 1, k = 1)
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1 - a0, a0, 0), nrow = n_mesh, ncol = 3, byrow = TRUE)
      pr0 <- pr0 * private$inside_
      pr0 <- pr0 / sum(pr0)
      return(pr0)
    },
    
    calc_D_llk = function() {
      D <- do.call(private$link2response_$D, list(self$par()$D))
      A <- sum(private$inside_) * private$dx_^2 / 1000^2
      n <- private$data_$n()
      pdet <- self$calc_pdet()
      llk <- n * log(D * A * pdet) - D * A * pdet - lfactorial(n)
      names(llk) <- NULL 
      return(llk)
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
      dt <- diff(self$data()$time())
      sd <- self$get_par("sd", m = 1)
      pdet <- C_calc_move_pdet(private$data_$n_occasions(), 
                               pr0, 
                               pr_empty, 
                               tpms, 
                               private$num_cells_,
                               private$inside_, 
                               private$dx_,
                               dt, 
                               sd,
                               3); 
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
      dt <- diff(self$data()$time())
      sd <- self$get_par("sd", m = 1)
      llk <- C_calc_move_llk(n, 
                             n_occasions,
                             pr0, 
                             pr_capture, 
                             tpms,
                             private$num_cells_, 
                             private$inside_, 
                             private$dx_, 
                             dt, 
                             sd, 
                             private$num_cores_, 
                             3, 
                             rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_pdet())
      llk <- llk + self$calc_D_llk()
      #plot(self$par()$beta[-1])
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    },
    
    simulate = function(seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      new_dat <- simulate_js_openscr(self$par(), 
                                     self$data()$n_occasions(), 
                                     self$data()$traps(), 
                                     self$data()$mesh(), 
                                     move = TRUE, 
                                     time = self$data()$time(), 
                                     seed = seed, 
                                     print = private$print_)
      return(new_dat)
    }
),
                   
  private = list(
    dx_ = NULL, 
    inside_ = NULL,
    num_cells_ = NULL,
    
    make_par = function() {
      samp_cov <- private$data_$covs(j = 1, k = 1, m = 1)
      n_par <- numeric(6)
      private$par_ <- vector(mode = "list", length = 6)
      for (par in 1:5) {
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
      private$par_[[6]] <- 0
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
        private$par_$sd[1] <-do.call(private$response2link_$sd, 
                                     list(start$sd))
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
      par$sd <- vec[grep("sd", names)]
      names(par$sd) <- gsub("sd.", "", names(par$sd))
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



