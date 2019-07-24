#' Setup regression spline design matrices for openpopscr 
#' 
#' This function is adapted from the gam.setup function in the R package mgcv. 
#' See ?mgcv:::gam.setup for more details.
#' 
#' This is a reduced form of gam.setup that only computed the design matrix.
#'
#' @param formula see ?mgcv::gam.setup
#' @param data see ?mgcv::gam.setup
#' @param knots see ?mgcv::gam.setup
#' @param sp see ?mgcv::gam.setup
#' @param min.sp see ?mgcv::gam.setup
#' @param H see ?mgcv::gam.setup
#' @param absorb.cons see ?mgcv::gam.setup
#' @param sparse.cons see ?mgcv::gam.setup
#' @param select see ?mgcv::gam.setup
#' @param idLinksBases see ?mgcv::gam.setup
#' @param scale.penalty see ?mgcv::gam.setup
#' @param paraPen see ?mgcv::gam.setup
#' @param gamm.call see ?mgcv::gam.setup
#' @param drop.intercept see ?mgcv::gam.setup
#' @param diagonal.penalty see ?mgcv::gam.setup
#' @param apply.by see ?mgcv::gam.setup
#' @param list.call see ?mgcv::gam.setup
#' @param modCon see ?mgcv::gam.setup
#'
#' @return A design matrix for the given formula and data 
#' @export
openpopscrgam <- function (formula, 
                           data = stop("No data supplied to gam.setup"), 
                           knots = NULL, 
                           sp = NULL, 
                           min.sp = NULL, 
                           H = NULL, 
                           absorb.cons = TRUE, 
                           sparse.cons = 0, 
                           select = FALSE, 
                           idLinksBases = TRUE, 
                           scale.penalty = TRUE, 
                           paraPen = NULL, 
                           gamm.call = FALSE, 
                           drop.intercept = FALSE, 
                           diagonal.penalty = FALSE, 
                           apply.by = TRUE, 
                           list.call = FALSE, 
                           modCon = 0) {
  # split formula into parametric and smooth parts
  split <- interpret.gam(formula)
  pterms <- terms(split$pf)
  # get number of smooths 
  if (length(split$smooth.spec) == 0) {
    if (split$pfok == 0) 
      stop("You've got no model....")
    m <- 0
  } else m <- length(split$smooth.spec)
  # compute parametric model frame 
  mf <- model.frame(split$pf, data, drop.unused.levels = FALSE)
  # drop intercept, if needed 
  if (drop.intercept) attr(pterms, "intercept") <- 1
  # get parametric design matrix 
  X <- model.matrix(pterms, mf)
  if (drop.intercept) {
    xat <- attributes(X)
    ind <- xat$assign > 0
    X <- X[, ind, drop = FALSE]
    xat$assign <- xat$assign[ind]
    xat$dimnames[[2]] <- xat$dimnames[[2]][ind]
    xat$dim[2] <- xat$dim[2] - 1
    attributes(X) <- xat
  }
  # get nonparametric design matrices 
  rownames(X) <- NULL
  sm <- list()
  newm <- 0
  if (m > 0) {
    for (i in 1:m) {
       # evaluate smooth 
        sml <- smoothCon(split$smooth.spec[[i]], data, 
                         knots, absorb.cons, scale.penalty = scale.penalty, 
                         null.space.penalty = select, sparse.cons = sparse.cons, 
                         diagonal.penalty = diagonal.penalty, apply.by = apply.by, 
                         modCon = modCon)
      # count how many smooth components 
      for (j in 1:length(sml)) {
        newm <- newm + 1
        sm[[newm]] <- sml[[j]]
      }
    }
    # add nonparametric to design matrix
    m <- newm
    for (i in 1:m) {
      X <- cbind2(X, sm[[i]]$X)
    }
  }
  return(X)
}
