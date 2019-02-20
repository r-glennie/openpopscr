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

#' Cormack-Jolly-Seber with transience model class 
#' 
#' @description Cormack-Jolly-Seber with transience model: fits model, formats inference, and 
#' simulates from fitted model. This model inherits all the functions of the 
#' CjsModel class. Functions that are different or additional are documented here. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item scr_data: a ScrData object 
#'   \item start: a named list of starting values 
#'   \item num_cores (optional, default = 1): number of processors cores to use 
#'   in parallelised code 
#'   \item print: (defualt TRUE) if TRUE then useful output is printed
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item simulate(): simulate ScrData object from fitted model
#' }
#' 
CjsTransientModel <- R6Class("CjsTransientModel", 
                         inherit = CjsModel, 
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
			index <- 1:data$n()
			if (print) cat("Computing entry occasions for each individual.......")
			private$entry_ <- apply(data$capthist(), 1, function(x) {min(index[rowSums(x) > 0])}) + 1 
			if (print) cat("done\n")
			if (print) cat("Reading formulae.......")
      private$form_ <- form 
      par_names <- sapply(form, function(f){f[[2]]})
      private$form_[[1]]<- form[par_names == "lambda0"][[1]]
      private$form_[[2]] <- form[par_names == "sigma"][[1]]
      private$form_[[3]] <- form[par_names == "phi"][[1]]
      private$form_[[4]] <- form[par_names == "sd"][[1]]
      private$form_ <- lapply(private$form_, function(f) {delete.response(terms(f))})
      names(private$form_) <- c("lambda0", "sigma", "phi", "sd")
      private$make_par() 
      private$link2response_ <- list(lambda0 = "exp", 
                            sigma = "exp", 
                            phi = "plogis", 
                            sd = "exp")
      private$response2link_ <- list(lambda0 = "log", 
                                     sigma = "log", 
                                     phi = "qlogis", 
                                     sd = "log") 
      if (print) cat("done\n") 
      if (print) cat("Initilising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$num_cores_ = num_cores
      private$print_ = print
    },
    
    calc_initial_distribution = function() {
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1, 0), nrow = n_mesh, ncol = 2, byrow = TRUE)
      pr0[, 1] <- pr0[, 1] * private$inside_ 
      pr0 <- pr0 / n_mesh
      return(pr0)
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
                             2, 
                             private$entry_)
      cat("llk:", llk, "\n")
      return(llk)
    },
    
    simulate = function(N = NULL, seed = NULL) {
      if (!is.null(N)) N <- self$data()$n()
      new_dat <- simulate_cjs_openscr(self$par(), 
                                      N, 
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
      n_par <- numeric(4)
      private$par_ <- vector(mode = "list", length = 3)
      for (par in 1:4) {
        X <- model.matrix(private$form_[[par]], data = samp_cov)
        n_par[par] <- ncol(X)
        if (par %in% c(3) & "t" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("t", private$data_$n_occasions() - 1, sep ="")]
        } else if (par %in% c(3) & "primary" %in% all.vars(private$form_[[par]])) {
          n_par[par] <- ncol(X) - 1
          par_vec <- rep(0, n_par[par])
          names(par_vec) <- colnames(X)[colnames(X) != paste("primary", private$data_$n_occasions() - 1, sep ="")]
        } else {
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
        private$par_$lambda0[1] <- do.call(private$response2link_$lambda0, 
                                           list(start$lambda0))
        private$par_$sigma[1] <-do.call(private$response2link_$sigma, 
                                           list(start$sigma))
        private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
        private$par_$sd[1] <-do.call(private$response2link_$sd, 
                                     list(start$sd))
        return(invisible())
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
      par$sd <- vec[grep("sd", names)]
      names(par$sd) <- gsub("sd.", "", names(par$sd))
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



