# logLik function 

#' Compute Log-Likelihood of Jolly-Seber Model
#'
#' @param object JsModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.JsModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(object$par()))
  class(llk) <- "logLik"
  return(llk)
}
