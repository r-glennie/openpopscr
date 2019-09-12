# Copyright (c) 2019 Richard Glennie
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

#' State model class 
#' 
#' @description Model for states individuals can move between during survey. 
#' \itemize{
#'   \item data : ScrData object 
#'   \item names: state names 
#'   \item structure: structure of state space graph in form of a edge matrix with 1 for estimated edge 
#'                    NA for edge with fixed transition probability, 0 for no edge 
#'   \item start: list with trm or tpm for initial transition rate matrix or transition probability matrix and 
#'         delta for initial distribution 
#'   \item cov: data frame, one row per state, one column per variable 
#'   \item delta_fixed: TRUE if initial is to be fixed, else it is estimated 
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k at mesh point m, similar to ScrData$covs() function
#' }
#' 
StateModel <- R6Class("StateModel", 
  public = list(
      initialize = function(data, names, structure, start, cov = NULL, delta_fixed = NULL) {
        private$nstates_ = length(names)
        private$names_ = names
        private$data_ = data
        private$check_struct(structure)
        private$struct_ = structure 
        private$delta_fixed_ = delta_fixed
        if (is.null(delta_fixed)) private$delta_fixed_ <- rep(FALSE, private$nstates_)
        if (is.null(cov)) {
          private$groups_ = data.frame(state = factor(1:private$nstates_))
        } else {
          private$groups_ <- cov 
          private$groups_$state = factor(1:private$nstates_)
        }
        private$make_par()
        private$initialise_par(start)
        private$compute_par()
      }, 
      
      nstates = function() {return(private$nstates_)}, 
      names = function(s = NULL) {
        if (is.null(s)) return(private$names_)
        return(private$names_[s])
      },
      par = function() {return(private$par_)}, 
      struct = function() {return(private$struct_)}, 
      delta = function() {return(private$delta_)},
      groups = function() {return(private$groups_)},
      groupnms = function() {return(names(private$groups_))}, 
      ngroups = function(c = 1) {return(nlevels(private$groups_[,c]))}, 
      
      trm = function(k = 1) {
        if (!identical(private$par_, private$last_par_)) private$compute_par() 
        private$last_par_ <- private$par_ 
        nstates <- self$nstates()
        Q <- matrix(0, nr  = nstates, nc = nstates)
        par <- private$computed_par_ 
        parloc <- private$parloc_
        if (nrow(parloc) > 0) {
          for (i in 1:nrow(parloc)) {
            Q[parloc[i,1], parloc[i, 2]] <- par[i, k]
          }
        }
        diag(Q) <- -rowSums(Q)
        return(Q)
      }, 
      
      tpm = function(k = 1, dt = 1) {
        Q <- self$trm(k)
        P <- expm(Q * dt)
        return(P)
      }, 
      
      plot_trm = function(k = 1, ...) {
        trm <- self$trm(k) 
        colnames(trm) <- self$names() 
        rownames(trm) <- self$names() 
        chain <- new("ctmc", generator = trm)
        plot(chain, ..., vertex.size = 30, edge.lty = 2)
        invisible()
      }, 
      
      plot_tpm = function(k = 1, ...) {
        tpm <- self$tpm(k) 
        colnames(tpm) <- self$names() 
        rownames(tpm) <- self$names() 
        chain <- new("markovchain", states = self$names(), transitionMatrix = tpm)
        plot(chain, ..., vertex.size = 30, edge.lty = 2)
        invisible()
      },
      
      set_par = function(par) {
        private$par_ <- par 
        private$compute_par()
      }, 
      
      estimates = function() {
        return(private$results_) 
      }, 
      
      print = function() {
        options(scipen = 999)
        if (is.null(private$results_)) {
          print("Attach to a model and then fit it.")
        } else {
          cat("STATE PARAMETER ESTIMATES (link scale)\n")
          print(signif(private$results_, 4))
        }
        options(scipen = 0)
      }, 
      
      calc_confint = function(est, sd, lcl, ucl) {
        private$results_ <- data.frame(est = est, sd = sd, lcl = lcl, ucl = ucl)
      }
      
  ), 
  private = list(
    nstates_ = NULL, 
    data_ = NULL, 
    names_ = NULL, 
    struct_ = NULL,
    nzeros_ = NULL, 
    parloc_ = NULL,
    parvecloc_ = NULL, 
    par_ = NULL,
    last_par_ = NULL, 
    computed_par_ = NULL, 
    Xmats_ = NULL,
    delta_ = NULL, 
    delta_fixed_ = NULL, 
    groups_ = NULL, 
    results_ = NULL, 
    
    check_struct = function(struct) {
      if (nrow(struct) != ncol(struct)) stop("Structure is not a square matrix.")
      if (any(diag(struct) != ".")) stop("Diagonal entries must be denoted with a full stop. They cannot 
                                         depend on covariates or be fixed as zero.")
      if (nrow(struct) != private$nstates_) stop("Dimensions of structure do not match length of names.")
      private$nzeros_ <- rowSums(struct == "0")
      if (any(private$nzeros_ > ncol(struct) - 1)) stop("Each row must have at least one non-zero entry.")
      cond <- (struct == "0") + (substr(struct, 0, 1) == "~") + (struct == ".")
      if (any(cond != 1)) stop("Structure entries must be characters with formulae starting with ~ or a zero or a full stop.")
      invisible()
    }, 
    
    make_par = function() {
      # covariates
      dat <- private$data_
      if (dat$n_primary() > 1) {
        covs <- as.data.frame(dat$covs(m = 1, j = 1, k = 1))
      } else {
        covs <- as.data.frame(dat$covs(m = 1, j = 1, p = 1))
      }
      struct <- self$struct() 
      nstates <- self$nstates() 
      nzeros <- private$nzeros_
      # get location of parameters in trm 
      parloc <- which(substr(struct, 0, 1) == "~", arr.ind = TRUE)
      # make X matrices and par matrix 
      npar <- rep(0, nrow(parloc))
      row <- NULL
      col <- NULL
      Xmats <- vector(mode = "list", length = nrow(parloc))
      if (nrow(parloc) > 0) {
        for (i in 1:nrow(parloc)) {
          Xmats[[i]] <- openpopscr:::openpopscrgam(as.formula(struct[parloc[i, 1], parloc[i, 2]]), data = covs)
          npar[i] <- ncol(Xmats[[i]])
          row <- c(row, rep(parloc[i,1], npar[i]))
          col <- c(col, rep(parloc[i,2], npar[i]))
        }
      }
      # add delta in 
      nz <- sum(private$delta_fixed_)
      ndeltapar <- nstates - nz - 1 
      if (ndeltapar > 0) {
        private$par_ <- rep(0, sum(npar) + ndeltapar)
        row <- c(rep(0, ndeltapar), row)
        col <- c(rep(0, ndeltapar), col)
      } else {
        private$par_ <- rep(0, sum(npar))
        ndeltapar <- 0 
      }
      private$parloc_ <- parloc 
      private$parvecloc_ <- matrix(0, nr = sum(npar) + ndeltapar, nc = 2)
      private$parvecloc_[,1] <- row 
      private$parvecloc_[,2] <- col
      private$Xmats_ <- Xmats 
      invisible()
    }, 
    
    initialise_par = function(start) {
      ## delta 
      if (is.null(start$delta)) stop("start$delta missing. Must supply starting values for delta.")
      nstates <- self$nstates()
      if (length(start$delta) < nstates) stop("start$delta must be of length == number of states")
      private$delta_ <- start$delta 
      nz <- sum(private$delta_fixed_)
      ndeltapar <- nstates - nz - 1 
      if (ndeltapar > 0) {
        private$par_[1:ndeltapar] <- invmlogit(start$delta[!private$delta_fixed_])
      }
      ## trm 
      if (!is.null(start$tpm) & is.null(start$trm)) {
        trm <- logm(start$tpm)
      } else {
        trm <- start$trm 
      }
      if (is.null(trm)) stop("Either start$trm or start$tpm must be given in starting values.")
      parloc <- private$parloc_ 
      parvecloc <- private$parvecloc_ 
      if (nrow(parloc) > 0) {
        for (i in 1:nrow(parloc)) {
          first <- min(which(parvecloc[,1] == parloc[i, 1] & parvecloc[,2] == parloc[i, 2]))
          private$par_[first] <- log(trm[parloc[i, 1], parloc[i, 2]])
        }
      }
    }, 
    
    compute_par = function() {
      parloc <- private$parloc_
      parvecloc <- private$parvecloc_ 
      par <- private$par_ 
      Xmats <- private$Xmats_
      pred <- matrix(0, nr = nrow(parloc), nc = private$data_$n_occasions())
      if (nrow(parloc) > 0) {
        for (i in 1:nrow(parloc)) {
          pred[i,] <- Xmats[[i]] %*% par[parvecloc[,1] == parloc[i,1] & parvecloc[,2] == parloc[i,2]]
        }
      }
      private$computed_par_ <- exp(pred) 
      nz <- sum(private$delta_fixed_)
      ndeltapar <- private$nstates_ - nz - 1 
      delta <- rep(0, private$nstates_)
      if (ndeltapar > 0) {
        delta[!private$delta_fixed_] <- mlogit(private$par_[1:ndeltapar])
        private$delta_ <- delta 
      } 
      invisible()
    }
  )
)
                      
                      
                      