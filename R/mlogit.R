#' Multinomial logit function 
#'
#' @param x vector to compute mlogit of or matrix where mlogit is applied by column  
#'
#' @return mlogit of vector x or each column of matrix x
#' @export
mlogit <- function(x) {
  y <- exp(x)
  mlogitVec <- function(v) {
    v <- c(1, v)
    v <- v / sum(v)
    return(v)
  }
  if (is.matrix(y)) {
    res <- apply(y, 2, mlogitVec)
  } else if (is.vector(y)) {
    res <- mlogitVec(y)
  } else if (is.array(y)) {
    if (length(dim(y)) < 2) stop("Arrays given to mlogit must have 2 or more dimensions.")
    res <- apply(y, 2:length(dim(y)), mlogitVec)
  } else {
    stop("Unrecognised object given to mlogit")
  }
  return(res)
}

#' Inverse Multinomial logit function
#'
#' @param y vector to compute inverse mlogit of or matrix where mlogit is applied by column  
#'
#' @return inverse mlogit of vector y or each column of matrix y 
#' @export
invmlogit <- function(y) {
  invmlogitVec <- function(v) {
    return(log(v[-1] / v[1]))
  }
  if (is.matrix(y)) {
    res <- apply(y, 2, invmlogitVec)
  } else if (is.vector(y)) {
    res <- invmlogitVec(y)
  } else if (is.array(y)) {
    if (length(dim(y)) < 2) stop("Arrays given to mlogit must have 2 or more dimensions.")
    res <- apply(y, 2:length(dim(y)), invmlogitVec)
  } else {
    stop("Unrecognised object given to invmlogit")
  }
  return(res)
}

