#' Temporal point process link function 
#'
#' @param x vector to compute pplink of or matrix where pplink is applied by column  
#'
#' @return pplink of vector x or each column of matrix x
#' @export
pplink <- function(x, dt) {
  y <- exp(x) * dt 
  pplinkVec <- function(v) {
    p <- v / sum(v)
    p1 <- exp(-sum(v))
    p <- (1 - p1) * p 
    p <- c(p1, p)
    return(p)
  }
  if (is.matrix(y)) {
    res <- apply(y, 2, pplinkVec)
  } else if (is.vector(y)) {
    res <- pplinkVec(y)
  } else if (is.array(y)) {
    if (length(dim(y)) < 2) stop("Arrays given to pplink must have 2 or more dimensions.")
    res <- apply(y, 2:length(dim(y)), pplinkVec)
  } else {
    stop("Unrecognised object given to pplink")
  }
  return(res)
}

#' Inverse temporal point process link function
#'
#' @param y vector to compute inverse pplink of or matrix where pplink is applied by column  
#'
#' @return inverse pplink of vector y or each column of matrix y 
#' @export
invpplink <- function(y, dt) {
  invpplinkVec <- function(v) {
    s <- -log(v[1])
    p <- v[-1] / (1 - v[1])
    p <- p * s 
    r <- log(p / dt)
    return(r)
  }
  if (is.matrix(y)) {
    res <- apply(y, 2, invpplinkVec)
  } else if (is.vector(y)) {
    res <- invpplinkVec(y)
  } else if (is.array(y)) {
    if (length(dim(y)) < 2) stop("Arrays given to pplink must have 2 or more dimensions.")
    res <- apply(y, 2:length(dim(y)), invpplinkVec)
  } else {
    stop("Unrecognised object given to invpplink")
  }
  return(res)
}

