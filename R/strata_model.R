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

#' Stratum model class
#' 
#' @description Fits the a model to different strata, sharing some parameters 
#' across strata. 
#' \itemize{
#'   \item scr_data: a list of ScrData objects, one per stratum  
#'   \item model: name of model class, e.g., ScrModel, CjsModel, JsModel
#'   \item shared_form: a named list of formulae for each parameter that is 
#'     shared across strata (~1 for constant)
#'   \item private_form: a named list of formula for each parameter that is private
#'   to each stratum, one per stratum
#'   \item print: (defualt TRUE) if TRUE then useful output is printed
#'   \item args : list of other arguments to passed to $new for model object 
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(...): calls get_par for each model object using arguments ... and
#'  returns list of results 
#'  \item calc_llk(): compute joint log-likelihood at current parameter values 
#'  \item fit: fit the model by obtaining the maximum likelihood estimates 
#'  \item par(): return current parameter of the model 
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates
#'  \item is_shared(par): returns indices of any named elements in par that 
#'  are shared parameters across strata 
#'  \item is_private(par): returns indices of any named elements in par that are 
#'  private parameters across strata  
#' }
#' 
StrataModel <- R6Class("StrataModel", 
  public = list(
    
    initialize = function(data, model, shared_form, private_form, start, print = TRUE, args = NULL) {
      private$data_ <- data
			private$n_strata_ <- length(data)
			args <- c(args, print = print)
			if (print) cat("Creating model objects for each stratum.........\n")
			private$objs_ <- private$forms_ <- vector(mode = "list", length = private$n_strata_)
			private$model_ <- model
			create_obj <- get(model)
			for (s in 1:private$n_strata_) {
			 if (print) cat("STRATUM", s, "\n")
			 private$forms_[[s]] <- c(shared_form, private_form[[s]]) 
			 private$objs_[[s]] <- do.call(create_obj$new, c(list(form = private$forms_[[s]], 
			                                 data = data[[s]], 
			                                 start = start), 
			                            args))
			}
			if (print) cat("Creating parameters for each model.......")
			private$shared_ <- sapply(shared_form, FUN = function(x){x[[2]]})
			private$private_ <- lapply(private_form, FUN = function(x){
			  sapply(x, FUN = function(y){y[[2]]})
			 })
			wpar <- private$convert_par2vec(private$objs_[[1]]$par())
	    private$par_ <- wpar[self$is_shared(wpar)]
	    private$npriv_ <- NULL
	    for (s in 1:private$n_strata_) {
	      wpar <- private$convert_par2vec(private$objs_[[s]]$par())
	      private$npriv_ <- c(private$npriv_, length(private$par_) + 1)
	      private$par_ <- c(private$par_, wpar[self$is_private(wpar)])
	    }
	    private$npriv_ <- c(private$npriv_, length(private$par_) + 1)
			if (print) cat("done\n")
      private$print_ = print
    },
    
    calc_llk = function(param = NULL, names = NULL) {
      if (!is.null(names)) names(param) <- names 
      if (!is.null(param)) {private$par_ <- param; private$split_par()}
      llk <- 0 
      for (i in 1:private$n_strata_) {
        ipar <- private$ipar_[[i]]
        nms <- names(ipar)
        llk <- llk + private$objs_[[i]]$calc_llk(ipar, nms)
      }
      return(llk)
    }, 
    
    get_object = function(i = NULL) {
      if (is.null(i)) return(private$objs_)
      return(private$objs_[[i]])
    },
    
    par = function() {
      return(private$par_)
    },
    
    get_par = function(...) {
      res <-  lapply(self$get_object(), FUN = function(m) {m$get_par(...)})
      return(res)
    },
    
    mle_llk = function() {
      res <- sapply(self$get_object(), FUN = AIC)
      return(sum(res))
    }, 
    
    fit = function() {
      par <- self$par()
      if (private$print_) cat("Fitting model..........\n")
      t0 <- Sys.time()
      mod <- suppressWarnings(optim(par, 
                                    private$calc_negllk,
                                    names = names(par),
                                    hessian = TRUE))
      t1 <- Sys.time()
      difft <- t1 - t0
      if (private$print_) cat("Completed model fitting in", difft, attr(difft, "units"), "\n")
      mle <- mod$par
      names(mle) <- names(par)
      code <- mod$convergence 
      if (code > 0) warning("model failed to converge with optim code ", code, "\n")
      if (private$print_ & code == 0) cat("Checking convergence.......converged", code, "\n")
      private$V_ <- solve(mod$hessian)
      private$llk_ <- -mod$value
      private$mle_ <- mle
      private$par_ <- mle
      private$split_par()
      private$split_V()
      private$set_mles()
    }, 
    
    is_shared = function(w) {
      names <- names(w)
      return(grep(paste(private$shared_, collapse = "|"), names))
    },
    
    is_private = function(w) {
      names <- names(w)
      ind <- grep(paste(private$shared_, collapse = "|"), names)
      return((1:length(w))[-ind])
    },
    
    print = function() {
      for (s in 1:private$n_strata_) {
        cat("Stratum ", s, "\n")
        print(private$objs_[[s]])
        cat("---------")
        cat("\n")
      }
    }
),
                   
  private = list(
    data_ = NULL,
    model_ = NULL, 
    n_strata_ = NULL, 
    forms_ = NULL,
    objs_ = NULL,
    shared_ = NULL, 
    private_ = NULL,
    par_ = NULL,
    ipar_ = NULL, 
    npriv_ = NULL,
    V_ = NULL, 
    iV_ = NULL, 
    llk_ = NULL,
    mle_ = NULL, 
		print_ = NULL, 
    
		calc_negllk = function(param = NULL, names = NULL) {
		  return(-self$calc_llk(param, names))
		}, 
		
		split_par = function() {
		  np <- private$npriv_
		  # shared parameters 
		  spar <- private$par_[1:(np[1] - 1)]
		  private$ipar_ <- vector(mode = "list", length = private$n_strata_)
		  for (s in 1:private$n_strata_) {
		    private$ipar_[[s]] <- c(spar, private$par_[np[s]:(np[s + 1] - 1)])
		  }
		  return(private$ipar_)
		}, 
		
		split_V = function() {
		  np <- private$npriv_
		  private$iV_ <- vector(mode = "list", length = private$n_strata_)
		  for (s in 1:private$n_strata_) {
		    private$iV_[[s]] <- private$V_[c(1:(np[1] - 1), np[s]:(np[s + 1] - 1)), 
		                                   c(1:(np[1] - 1), np[s]:(np[s + 1] - 1))]
		  }
		  return(private$ipar_)
		}, 
		
		set_mles = function() {
		  for (s in 1:private$n_strata_) {
		    private$objs_[[s]]$set_mle(private$ipar_[[s]], private$iV_[[s]], private$llk_)
		  }
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
		}
  )                 
)



