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
#'   \item print (default = TRUE): if TRUE then helpful output is printed to the screen
#' }
#' 
#' Methods include those in JsModel with the following overwritten: 
#' \itemize{
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#' }
#' 
JsTransientModel <- R6Class("JsTransientModel", 
                            inherit = JsModel, 
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      if (print) cat("Creating rectangular mesh......")
      newmesh <- rectangularMask(data$mesh())
      private$dx_ <- attr(newmesh, "spacing")
      private$inside_ <- as.numeric(pointsInPolygon(newmesh, data$mesh()))
      cov_list <- data$get_cov_list() 
      private$data_ <- ScrData$new(data$capthist(), newmesh, data$time(), data$primary())
      covnms <- names(cov_list$cov)
      allcovs <- names(private$data_$get_cov_list()$cov)
      covnms <- covnms[!(covnms %in% allcovs)]
      if (!is.null(covnms)) {
        for (i in seq(covnms)) {
          private$data_$add_covariate(covnms[i], cov_list$cov[[i]], cov_list$cov_type[i])
        }
      }
      box <- attributes(newmesh)$boundingbox
      region <- c(diff(box[1:2, 1]), diff(box[c(1, 3), 2]))
      private$num_cells_ <- numeric(3)
      private$num_cells_[1] <- nrow(newmesh)
      private$num_cells_[2] <- floor(region[1] / private$dx_)
      private$num_cells_[3] <- nrow(newmesh) / private$num_cells_[2]
      if (print) cat("done\n")
      if (print) cat("Reading formulae.......")
      order <- c("phi", "beta", "sd", "D")
      private$read_formula(form, detectfn, order)
      # add other parameters type
      private$par_type_[private$detfn_$npars() + 1] <- "km"
      private$par_type_[private$detfn_$npars() + 2] <- "k1m"
      private$par_type_[private$detfn_$npars() + 3] <- "km"
      private$par_type_[private$detfn_$npars() + 4] <- "m"
      names(private$form_) <- c(private$detfn_$pars(), "phi", "beta", "sd", "D")
      # make parameter list 
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"), list("mlogit"), list("exp"), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi", "beta", "sd", "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"), list("invmlogit"), list("log"), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi", "beta", "sd", "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      if (print) cat("done\n")
      private$print_ = print 
    },
    
    calc_initial_distribution = function() {
      a0 <- self$get_par("beta", j = 1, k = 1)
      n_mesh <- private$data_$n_meshpts()
      pr0 <- matrix(c(1 - a0, a0, 0), nrow = n_mesh, ncol = 3, byrow = TRUE)
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:n_mesh) * a * private$inside_
      pr0[,1] <- pr0[,1] * D
      pr0[,2] <- pr0[,2] * D 
      return(pr0)
    },
    
    calc_D_llk = function() {
      n <- private$data_$n()
      Dpdet <- self$calc_Dpdet()
      llk <- n * log(sum(Dpdet)) - sum(Dpdet) - lfactorial(n)
      names(llk) <- NULL 
      return(llk)
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
      dt <- diff(self$data()$time())
      sd <- self$get_par("sd", m = 1)
      Dpdet <- C_calc_move_pdet(private$data_$n_occasions(), 
                               pr0, 
                               pr_empty, 
                               tpms, 
                               private$num_cells_,
                               private$inside_, 
                               private$dx_,
                               dt, 
                               sd,
                               3); 
      a <- private$data_$cell_area()
      D <- self$get_par("D", m = 1:private$data_$n_meshpts()) * a * private$inside_
      Dpdet <- sum(D) - Dpdet
      return(Dpdet)
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
                             3, 
                             rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_Dpdet())
      llk <- llk + self$calc_D_llk()
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    }
),
                   
  private = list(
    dx_ = NULL, 
    inside_ = NULL,
    num_cells_ = NULL,
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                               list(start[[names[i]]]))
      }
      private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
      private$par_$beta[1] <- do.call(private$response2link_$beta,
                                      list(c(start$beta, rep((1 - start$beta) / (self$data()$n_occasions() - 1), self$data()$n_occasions() - 1))))[1]
      private$par_$sd[1] <-do.call(private$response2link_$sd, 
                                     list(start$sd))
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                         list(start$D / private$data_$n_meshpts()))
      # compute initial parameters for each jkm
      private$compute_par()
      return(invisible())
    }
  )                 
)



