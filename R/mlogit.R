#' Multinomial logit function 
#'
#' @param x vector to compute mlogit of or matrix where mlogit is applied by column  
#'
#' @return mlogit of vector x or each column of vector x
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
  } else {
    res <- mlogitVec(y)
  }
  return(res)
}

#' Inverse Multinomial logit function
#'
#' @param y vector to compute inverse mlogit of or matrix where mlogit is applied by column  
#'
#' @return inverse mlogit of vector y or each column of vector y 
#' @export
invmlogit <- function(y) {
  invmlogitVec <- function(v) {
    return(log(v[-1] / v[1]))
  }
  if (is.matrix(y)) {
    res <- apply(y, 2, invmlogitVec)
  } else {
    res <- invmlogitVec(y)
  }
  return(res)
}

