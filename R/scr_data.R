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
#'     (used for irregularly spaced occasions). Default is vector 1:n_occasions.
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
#'  \item plot_mesh(cov, nbreaks = 10, ...): if cov not NULL then plots covariate value over mesh, nbreaks 
#'        is number of breaks in covariate value plotted, ... is passed to geom_point calls
#'  \item detector_type(): return detector type index (1 = count, 2 = proximity, 
#'        3 = multi/single)
#'  \item get_cov_list(): return list of covariates
#'  \item covs(j, k, m): return covariate values at detector j 
#'        on occasion k at mesh point m, if any of j,k are NULL then returns for all values; if m is NULL returns spatial mean 
#'  \item n(): return number of individuals seen over the entire survey
#'  \item n_occasions(): return number of primary occasions in the survey, n_occasions("all") returns #'        total number of secondary occasions, n_occasions(i) returns number of secondary occasions #'        in primary occasion i 
#'  \item n_traps(): return number of traps used at some time in the survey
#'  \item n_meshpts(): return number of mesh points 
#'  \item n_primary(): return number of primary occasions 
#'  \item n_secondary(): return vector of number of secondary occasions per primary
#'  \item primary(): return vector showing primary occasion each secondary occasion is assigned to
#'  \item encrate(each = FALSE): return empirical estimate of encounter rate, each = TRUE returns for each primary occasion
#'  \item encrange(k = NULL, each = FALSE): return empirical estimate of encounter range, each = TRUE returns for each primary occasion, k can be a vector/numeric specifying which occasions to compute over
#'  \item unique: return number of unique captures per occasion 
#'  \item area(): return total area of the mesh 
#'  \item cell_area(): return area of a single grid cell on mesh
#'  \item replace_mesh(newmesh): replace stored mesh with newmesh
#'  \item distances(): return matrix with (i,j)th entry the distance from trap i to mesh point j
#'  \item add_covariate(covariate_name, covariate_vector/matrix, covariate_type): add covariate to object
#'  \item remove_covariate(covariate_name): remove covariate from data object
#' }
#' 
ScrData <- R6Class("ScrData", 
  public = list( 
    initialize = function(capthist, 
                          mesh, 
                          time = NULL, 
                          primary = NULL) {
      private$check_input(capthist, mesh, time, primary) 
      ## detectors
      private$detector_type_ <- switch(attr(traps(capthist), "detector")[1], 
                                       count = 1, 
                                       proximity = 2, 
                                       multi = 3,
                                       single = 3, 
                                       transect = 4, 
                                       transectX = 5, 
                                       polygon = 6, 
                                       polygonX = 7)
      if (attr(traps(capthist), "detector")[1] == "single") warning("Single-catch detectors treated as multi-catch detectors.")
      private$capthist_ <- capthist 
      if (is.null(secr::usage(secr::traps(capthist)))) {
        if (private$detector_type_ %in% 4:7) {
          nids <- length(unique(attr(traps(capthist), "polyID")))
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = nids, nc = dim(capthist)[2])
        } else {
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = dim(capthist)[3], nc = dim(capthist)[2])
        }
      }
      if (private$detector_type_ %in% 4:7) {
        outputdetector <- ifelse(private$detector_type_ %in% c(4, 6), "count", "proximity")
        spacing <- attr(mesh, "polygon_trap_spacing")
        if (is.null(spacing)) spacing <- attr(mesh, "spacing")
        private$capthist_ <- secr::discretize(capthist, 
                                              spacing = spacing, 
                                              outputdetector = outputdetector, 
                                              cell.overlap = TRUE)
        if (outputdetector == "proximity") attributes(traps(private$capthist_))$detector <- "multi"
        private$detector_type_ <- ifelse(outputdetector == "count", 1, 3)
      }
      ## split capthist into primary occasions 
      if (is.null(primary)) primary <- rep(1, dim(capthist)[2])
      private$primary_ <- primary 
      private$n_primary_ <- max(primary)
      private$capthists_ <- secr:::split.capthist(capthist, as.factor(primary), byoccasion = TRUE, dropnullCH = FALSE, dropunused = FALSE)
      private$n_occasions_ <- length(private$capthists_)
      if (private$n_occasions_ == 1) private$n_occasions_ <- dim(capthist)[2]
      private$mesh_ <- mesh
      if (is.null(time)) {
        private$time_ <- seq(1, self$n_occasions("all"))
      } else {
        private$time_ <- time
      }
      ## add built-in covariates
      private$cov_$t <- (1:self$n_occasions("all") - 1)
      private$cov_$T <- as.factor((1:self$n_occasions("all") - 1))
      private$cov_type_ <- c("k", "k")     
      private$cov_$primary <- (private$primary_ - 1)
      private$cov_$Primary <- as.factor((private$primary_ - 1))
      private$cov_type_ <- c(private$cov_type_, "k", "k")   
      private$cov_$x <- scale(private$mesh_[,1])[,1]
      private$cov_$y <- scale(private$mesh_[,2])[,1]
      private$cov_type_ <- c(private$cov_type_, "m")
      private$cov_type_ <- c(private$cov_type_, "m")
      if (!(private$detector_type_ %in% 1:7)) stop("Detector type not implemented.")
      ## compute distance trap-to-mesh
      self$calc_distances()
    },
    
    #### OUTPUT FUNCTIONS 
    
    print = function(i = ".") {
       plot(self$mesh())
       varycol <- ifelse(self$n() > 1000, FALSE, TRUE)
       plot(self$capthist(i), add = T, varycol = varycol)
       plot(self$traps(), add = T)
       print(summary(self$capthist())[[4]])
    },
    
    plot_mesh = function(cov = NULL, nbreaks = 10, ...) {
      mesh <- self$mesh()
      mesh[,1] <- scale(mesh[,1], scale = F)
      mesh[,2] <- scale(mesh[,2], scale = F)
      meshdat <- data.frame(x = mesh[,1], y = mesh[,2])
      var <- NULL
      if (!is.null(cov)) {
        if (is.character(cov)) {
          nmesh <- self$n_meshpts()
          if (cov %in% names(self$covs())) {
            f <- which(names(self$covs()) == cov) 
            if (!grepl("m", private$cov_type_[f])) stop ("Covariate must be spatial.")
            var <- self$covs()[cov][[1]]
            name <- cov 
          }
        } else {
          var <- cov 
          name <- ""
        }
        if (is.factor(var)) nbreaks <- nlevels(var)
        collist <- viridis(nbreaks)
        fvar <- cut(as.numeric(var), breaks = nbreaks)
        meshdat$cols <- collist[fvar] 
      }
      if (is.factor(var)) {
        plt <- ggplot(meshdat) + geom_point(aes(x = x, y = y, col = cols), ...) + 
          scale_color_viridis_d("Legend", labels = levels(var)) 
      } else if (!is.null(var)) {
        plt <- ggplot(meshdat) + geom_point(aes(x = x, y = y, col = cols), ...) + 
          scale_color_viridis_d("Value", labels = levels(fvar)) 
      } else {
        plt <- ggplot(meshdat) + geom_point(aes(x = x, y = y), col = "grey", ...) +
          theme(legend.position = "none")
      }
      plt <- plt + theme_bw()
      return(plt)
    }, 

    #### ACCESSORS 
      
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
    
    covs = function(j = NULL, k = NULL, m = NULL, i = NULL) {
      if (any(is.null(j))) j0 <- seq(1, self$n_traps()) else j0 <- j
      if (any(is.null(k))) k0 <- seq(1, self$n_occasions("all")) else k0 <- k 
      if (any(is.null(m))) m0 <- seq(1, self$n_meshpts()) else m0 <- m 
      if (any(is.null(i))) i0 <- 1 else i0 <- i 
      if (!is.null(i) && length(i) > 1) stop("For covs(), arguument i must be a single integer.")
      if (!is.null(i) && length(k) > 1) stop("For covs(), argument k must be a single integer when i is used.")
      dat <- lapply(1:length(private$cov_), FUN =  function(c) {
        switch(private$cov_type_[c], 
               j = private$cov_[[c]][j0], 
               k = private$cov_[[c]][k0], 
               m = private$cov_[[c]][m0], 
               kj = private$cov_[[c]][k0, j0],
               km = private$cov_[[c]][k0, m0], 
               i = private$cov_[[c]][[i0]], 
               ik = private$cov__[[c]][[i0]][[k0]], 
               private$cov_[[c]])
      })
      names(dat) <- names(private$cov_)
      return(dat)
    }, 
    
    n = function(i = ".") {
      if (i == ".") return(dim(private$capthist_)[1])
      return(dim(private$capthists_[[i]]))
    },
    n_occasions = function(i = NULL) {
      if (!is.null(i)) {
        if (i == "all") {
          return(dim(private$capthist_)[2])
        } else {
          return(sum(self$primary() == i))
        }
      } else {
        return(private$n_occasions_)
      }
    },
    n_traps = function() {
      # number of transects
      if (private$detector_type_ %in% 4:7) {
        return(length(unique(attr(traps(private$capthist_),"polyID"))))
      }
      return(dim(private$capthist_)[3])
    }, 
    n_meshpts = function() {return(dim(private$mesh_)[1])},
    n_primary = function() {return(private$n_primary_)}, 
    n_secondary = function() {return(as.numeric(table(self$primary())))}, 
    primary = function() {return(private$primary_)},
    area = function() {return(self$n_meshpts() * attributes(private$mesh_)$area * 0.01)},
    cell_area = function() {return(attributes(private$mesh_)$area * 0.01)}, 
    distances = function(){return(private$distances_)},
    
    #### SUMMARY STATISTICS 
    
    encrate = function(each = FALSE) {
      ndet <- as.numeric(summary(self$capthist())[[4]][6,1:self$n_occasions("all")])
      nind <- as.numeric(summary(self$capthist())[[4]][1,1:self$n_occasions("all")])
      enc <- ndet/nind
      if (each) {
        return(enc)
      } else {
        return(mean(enc))
      }
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
    unique = function() {
      return(as.numeric(summary(self$capthist())[[4]][2,1:self$n_occasions("all")]))
    },
    
    #### FUNCTIONS 
    replace_mesh = function(newmesh) {
      if (!("mask" %in% class(newmesh))) stop("Invalid mesh object.")
      private$mesh_ <- newmesh 
      self$calc_distances()
    }, 
    
    calc_distances = function() {
      private$distances_ <- t(apply(self$traps(), 1, private$dist_to_row))
    }, 
        
    add_covariate = function(cov_name, cov, cov_type) {
      names <- names(private$cov_)
      ncov <- length(private$cov_)
      if (!is.character(cov_name)) stop("Covariate name must be a character string.")
      if (cov_name %in% names) stop("Covariate with that name already exists.")
      if (!(cov_type %in% c("i", "ik", "j", "k", "kj", "km", "m"))) stop("Invalid covariate type.")
      if (!(cov_type %in% c("i", "ik"))) {
        if (!is.factor(cov) & !is.numeric(cov)) stop("Invalid covariate, must be factor or numeric.")
      }
      if (cov_type == "i") if (length(cov) != self$n()) stop("Invalid covariate.")
      if (cov_type == "ik") if (nrow(cov) != self$n() || ncol(cov) != self$n_occasions()) stop("Invalid covariate.")
      if ("i" %in% private$cov_type_ || "ik" %in% private$cov_type_) stop("Can only have one known state covariate.")
      if (cov_type == "j") if (length(cov) != self$n_traps()) stop("Invalid covariate.")
      if (cov_type == "k") if (length(cov) != self$n_occasions()) stop("Invalid covariate.")
      if (cov_type == "kj") if (nrow(cov) != self$n_occasions() || ncol(cov) != self$n_traps()) stop("Invalid covariate.")
      if (cov_type == "km") if (nrow(cov) != self$n_occasions() || ncol(cov) != self$n_meshpts()) stop("Invalid covariate.")
      if (cov_type == "m") if (length(cov) != self$n_meshpts()) stop("Invalid covariate.")
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
    capthist_ = NULL, # capture history object 
    capthists_ = NULL, # list of capthists split into primary occasions 
    mesh_ = NULL, # mesh object 
    time_ = NULL, # vector of occasion start times 
    cov_ = NULL, # list of covariates 
    cov_type_ = NULL, # type of each covariate in cov_ list 
    distances_ = NULL, # matrix of distances from trap to mesh 
    detector_type_ = NULL, # type of detectors as integer (see initialize)
    primary_ = NULL, # vector of indices for each occasion indexing what primary it belongs to
    n_primary_ = NULL, # number of primary occasions
    n_occasions_ = NULL, # number of occasions
    
    #### FUNCTIONS
    
    # used to compute distances from trap-to-mesh 
    dist_to_row = function(r) {
      dist <- function(r2) {
        sqrt(sum((r2 - r)^2))
      }
      apply(private$mesh_, 1, dist)
    }, 
    
    # check input into intialize 
    check_input = function(capthist, mesh, time, primary) {
      if (!("capthist" %in% class(capthist))) stop("Invalid capture history object.")
      if (!("mask" %in% class(mesh))) stop("Invalid mesh object.")
      if (!is.null(time)) {
        if (length(time) != dim(capthist)[2]) stop("Length of time vector not equal to number of occasions.")
        if (!is.numeric(time) | !is.vector(time)) stop("Time is not a numeric vector.")
        if (max(abs(sort(time) - time)) > 1e-10) stop("Time must be an increasing vector of numbers.")
      }
      if (!is.null(primary)) {
        if (!is.numeric(primary) | !is.vector(primary)) stop("Primary is not a numeric vector.")
        if (length(primary) != dim(capthist)[2]) stop("Length of primary vector not equal to number of occasions.")
        if (max(abs(round(primary,0) - primary)) > 1e-10) stop("Primary must be integer labels.")
        if (max(abs(sort(primary) - primary)) > 1e-10) stop("Primary must be an increasing vector of numbers.")
        nprim <- max(primary)
        testprim <- sort(unique(primary))
        if (length(testprim) != nprim) stop("Primary must contain integers 1:maximum primary.")
        if (max(abs(testprim - 1:nprim)) > 0.5) stop("Primary must contain integers 1:maximum primary.")
      }
      return(0)
    }
  )
)
