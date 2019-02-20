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
# file description: scr data class based on secr package objects 
#
################################################################################

#' Spatial capture-recapture data class 
#' 
#' @description Encapsulates data from spatial capture-recapture survey into a single 
#' object. Object can be created using $new with arguments: 
#' \itemize{
#'   \item capthist: capture history object from secr package with trap information included
#'   \item mesh: mask object from secr package
#'   \item time: optional vector of numeric time each occasion took place at 
#'     (used for irregularly spaced occasions). Default is vector of ones.
#'   \item cov: a list of covariate (for detector and occasion)
#'   \item cov_type: type of covariate (j = detector cov, k = occasion cov, jk = detector+occasion,)
#'   \item primary: vector with index for each occasion in capthist that pools occasions into primary occasions 
#'        
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item capthist(): return capture history object 
#'  \item traps(): return traps object 
#'  \item mesh(): return mesh object
#'  \item time(): return time vector 
#'  \item covs(j, k): return covariate values at detector j 
#'        on occasion k  
#'  \item n(): return number of individuals seen over the entire survey
#'  \item n_occasions(): return number of capture occasions in the survey
#'  \item n_traps(): return number of traps used at some time in the survey
#'  \item n_meshpts(): return number of mesh points 
#'  \item area(): return total area of the mesh 
#'  \item add_covariate(covariate_name, covariate_dataframe, covariate_type): add covariate to object
#'  \item distance(): return matrix with (i,j)th entry the distance from trap i to mesh point j
#' }
#' 
ScrData <- R6Class("ScrData", 
  public = list( 
    initialize = function(capthist, 
                          mesh, 
                          time = NULL, 
                          cov = NULL, 
                          cov_type = NULL, 
                          primary = NULL) {
      private$capthist_ <- capthist 
      ## split capthist into primary occasions 
      if (is.null(primary)) primary <- rep(1, dim(capthist)[2])
      private$primary_ <- primary 
      private$n_primary_ <- max(primary)
      private$capthists_ <- secr:::split.capthist(capthist, as.factor(primary), byoccasion = TRUE, dropnullCH = FALSE, dropunused = FALSE)
      private$n_occasions_ <- length(private$capthists_)
      if (private$n_occasions_ == 1) private$n_occasions_ <- dim(capthist)[2]
      private$detector_type_ <- switch(attr(traps(capthist), "detector")[1], 
                                       count = 1, 
                                       proximity = 2, 
                                       multi = 3,
                                       transect = 4)
      if (is.null(secr::usage(secr::traps(capthist)))) {
        if (private$detector_type_ == 4) {
          nids <- length(unique(attr(traps(capthist), "polyID")))
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = nids, nc = dim(capthist)[2])
        } else {
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = dim(capthist)[3], nc = dim(capthist)[2])
        }
      }
      private$mesh_ <- mesh
      if (is.null(time)) {
        private$time_ <- seq(1, self$n_occasions("all"))
      } else {
        private$time_ <- time
      }
      private$cov_ <- cov 
      private$cov_$t <- as.factor((1:self$n_occasions("all") - 1))
      private$cov_type_ <- c(cov_type, "k")
      
      private$cov_$primary <- as.factor((private$primary_ - 1))
      private$cov_type_ <- c(private$cov_type_, "k")
      
      if (!(private$detector_type_ %in% c(1,2,3,4))) stop("openpopscr only implements 'count', 'proximity', 'multi', and 'transect' detectors, 
you are using detectors of another type.")
      self$calc_distances()
    },
    
    print = function(i = ".") {
       plot(self$mesh())
       varycol <- ifelse(self$n() > 1000, FALSE, TRUE)
       plot(self$capthist(i), add = T, varycol = varycol)
       plot(self$traps(), add = T)
       print(summary(self$capthist())[[4]])
    },
   
    capthist = function(i = ".") {
      if (i == ".") return(private$capthist_)
      return(private$capthists_[[i]])
    },
    traps = function(i = ".") {
      if (i == ".") return(traps(private$capthist_))
     return(traps(private$capthists_[[i]]))
    },
    mesh = function() {return(private$mesh_)},
    time = function() {return(private$time_)},
    detector_type = function() {return(private$detector_type_)}, 
    get_cov_list = function() {return(list(cov = private$cov_, cov_type = private$cov_type_))}, 
    
    covs = function(j = NULL, k = NULL, m = NULL) {
      if (any(is.null(j))) j0 <- seq(1, self$n_traps()) else j0 <- j
      if (any(is.null(k))) k0 <- seq(1, self$n_occasions("all")) else k0 <- k 
      if (any(is.null(m))) m0 <- seq(1, self$n_meshpts()) else m0 <- m 
      dat <- lapply(1:length(private$cov_), FUN =  function(c) {
        switch(private$cov_type_[c], 
               j = private$cov_[[c]][j0], 
               k = private$cov_[[c]][k0], 
               m = private$cov_[[c]][m0], 
               jk = private$cov_[[c]][j0, k0], 
               jm = private$cov_[[c]][j0, m0], 
               km = private$cov_[[c]][k0, m0], 
               jkm = private$cov_[[c]][j0, k0, m0], 
               private$cov_[[c]])
      })
      names(dat) <- names(private$cov_)
      return(dat)
    }, 
    
    n = function(i = ".") {
      if (i == ".") return(dim(private$capthist_)[1])
      return(dim(private$capthists_[[i]]))
    },
    n_occasions = function(i = ".") {
      if (i == ".") return(private$n_occasions_)
      if (i == "all") return(dim(private$capthist_)[2])
      return(dim(private$capthists_[[i]])[2])
    },
    n_traps = function() {
      # number of transects
      if (private$detector_type_ == 4) return(length(unique(attr(traps(private$capthist_), "polyID"))))
      return(dim(private$capthist_)[3])
      }, 
    n_meshpts = function() {return(dim(private$mesh_)[1])},
    n_primary = function() {return(private$n_primary_)}, 
    n_secondary = function() {return(as.numeric(table(self$primary())))}, 
    primary = function() {return(private$primary_)},
    encrate = function() {
      ndetections <- summary(self$capthist())[[4]][1, self$n_occasions(".") + 1]
      rate <- ndetections / (self$n() * self$n_occasions("."))
      return(rate)
    },
    encrange = function(k = NULL, each = FALSE) {
      if (is.null(k)) k <- seq(1, self$n_occasions())
      if (!each | length(k) == 1) {
        if (private$n_primary_ == 1) {
          subcap <- subset(self$capthist(), occasions = k)
        } else {
          wh <- which(self$primary() %in% k) 
          subcap <- subset(self$capthist(), occasions = wh)
        }
        return(RPSV(subcap))
      } else {
        if (private$n_primary_ == 1) {
          range <- numeric(length(k))
          for (i in seq(k)) range[i] <- RPSV(subset(self$capthist(), occasions = k[i]))
        } else {
          range <- numeric(length(k))
          for (i in seq(k)) range[i] <- RPSV(self$capthist(i))
        }
        return(range)
      }
    },
    
    area = function() {return(self$n_meshpts() * attributes(private$mesh_)$area * 0.01)},
    distances = function(){return(private$distances_)},
    
    calc_distances = function() {
      dist_to_row <- function(r) {
        dist <- function(r2) {
          sqrt(sum((r2 - r)^2))
        }
        apply(private$mesh_, 1, dist)
      }
      private$distances_ <- t(apply(self$traps(), 1, dist_to_row))
      ## if transect detectors, then take minimum distance approached
      if (self$detector_type() == 4) {
        ids <- attr(self$traps(), "polyID")
        unique_ids <- unique(ids)
        nids <- length(unique_ids)
        newdist <- matrix(0, nr = nids, nc = self$n_meshpts())
        for (i in 1:nids) {
          subdist <- private$distances_[ids == unique_ids[i],]
          newdist[i,] <- apply(subdist, 2, min)
        }
        private$distances_ <- newdist 
      } 
    }, 
    
    add_covariate = function(cov_name, cov, cov_type) {
      names <- names(private$cov_)
      ncov <- length(private$cov_)
      private$cov_[[ncov + 1]] <- cov 
      names(private$cov_) <- c(names, cov_name)
      private$cov_type_ <- c(private$cov_type_, cov_type)
    },
    
    remove_covariate = function(cov_name) {
      names <- names(self$covs())
      id <- which(names == cov_name)
      if (length(id) == 0) stop("No covariates with that name.")
      private$cov_ <- subset(private$cov_, names != cov_name)
      names(private$cov_) <- names[-id]
      private$cov_type_ <- private$cov_type_[-id]
    }
    
  ), 
  
  private = list(
    capthist_ = NULL, 
    capthists_ = NULL, 
    mesh_ = NULL,
    time_ = NULL,
    cov_ = NULL, 
    cov_type_ = NULL, 
    distances_ = NULL,
    detector_type_ = NULL, 
    primary_ = NULL, 
    n_primary_ = NULL, 
    n_occasions_ = NULL
  )
)








