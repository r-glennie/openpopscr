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
#'   \item print: (default TRUE) if TRUE then useful output is printed
#' }
#' 
#' Methods include those in CjsModel with the following overwritten: 
#' \itemize{
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#' }
#' 
CjsTransientModel <- R6Class("CjsTransientModel", 
                         inherit = CjsModel, 
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, statemod = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      private$data_ <- data
      private$start_ <- start 
      private$dx_ <- attr(data$mesh(), "spacing")
      private$inside_ <- matrix(-1, nr = data$n_meshpts(), nc = 4)
      for (m in 1:data$n_meshpts()) {
        dis <- sqrt((data$mesh()[m, 1] - data$mesh()[,1])^2 + (data$mesh()[m, 2] - data$mesh()[, 2])^2)
        wh <- which(dis < (1 + 1e-6) * private$dx_ & dis > 1e-16) - 1 
        if(length(wh) > 0) private$inside_[m, 1:length(wh)] <- as.numeric(wh) 
      }
      box <- attributes(data$mesh())$boundingbox
      region <- c(diff(box[1:2, 1]), diff(box[c(1, 3), 2]))
      private$num_cells_ <- numeric(3)
      private$num_cells_[1] <- data$n_meshpts()
      private$num_cells_[2] <- floor(region[1] / private$dx_)
      private$num_cells_[3] <- data$n_meshpts() / private$num_cells_[2]
			index <- 1:data$n()
			if (print) cat("Computing entry occasions for each individual.......")
			private$entry_ <- apply(data$capthist(), 1, function(x) {min(index[rowSums(x) > 0])}) - 1
			if (private$data_$n_primary() > 1) {
			  private$entry_ <- private$data_$primary()[private$entry_ + 1] - 1 
			}
			if (print) cat("done\n")
			if (print) cat("Reading formulae.......")
			order <- c("phi", "sd")
			private$read_formula(form, detectfn, statemod, order)
			# add parameters other than detection 
			private$par_type_[private$detfn_$npars() + 1] <- "p1ms"
			private$par_type_[private$detfn_$npars() + 2] <- "p1ms"
			names(private$form_) <- c(private$detfn_$pars(), "phi", "sd")
			# create parameter list
      private$make_par() 
      private$link2response_ <- c(private$detfn_$link2response(), list("plogis"), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "phi", "sd")
      private$response2link_ <- c(private$detfn_$response2link(), list("qlogis"), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "phi", "sd")
      if (print) cat("done\n") 
      if (print) cat("Initilising parameters.......")
      private$initialise_par(start)
      private$read_states() 
      if (print) cat("done\n")
      private$print_ = print
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
      # get tpms for state model 
      nstates <- self$state()$nstates() + 1 
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
      sd <- self$get_par("sd", s = 1:self$state()$nstates())
      sd[is.na(sd)] <- -10
      if (any(sd > 10 * private$start_$sd)) return(-Inf)
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
                             nstates,
                             0, 
                             1,
                             private$entry_)
      # compute probability of initial detection
      inipdet <- self$calc_initial_pdet(pr_capture) 
      llk <- llk - sum(log(inipdet))
      if(private$print_) cat("llk:", llk, "\n")
      return(llk)
    }
),
                   
  private = list(
		dx_ = NULL, 
		inside_ = NULL,
		num_cells_ = NULL,
		start_ = NULL, 
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                        list(start[[names[i]]]))
      }
      names(private$par_) <- c(names, "phi", "sd")
      private$par_$phi[1] <- do.call(private$response2link_$phi, 
                                           list(start$phi))
      private$par_$sd[1] <-do.call(private$response2link_$sd, 
                                     list(start$sd))
      # compute initial parameters for each jkm
      private$compute_par()
        return(invisible())
    }
  )                 
)



