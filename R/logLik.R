# logLik function 

#' Compute Log-Likelihood of SCR Model
#'
#' @param object ScrModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.ScrModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of SCR Transient Model
#'
#' @param object ScrTransientModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.ScrTransientModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of Cormack-Jolly-Seber Model
#'
#' @param object CjsModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.CjsModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of Cormack-Jolly-Seber Transient Model
#'
#' @param object CjsModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.CjsTransientModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of Jolly-Seber Model
#'
#' @param object JsModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.JsModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of Jolly-Seber Transient Model
#'
#' @param object JsModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.JsTransientModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}

#' Compute Log-Likelihood of Strata Model
#'
#' @param object StrataModel object 
#' @param ... ignored  
#'
#' @return log-likelihood 
#' @export
logLik.StrataModel <- function(object, ...) {
  llk <- object$mle_llk()
  attributes(llk) <- list(df = length(unlist(object$par())))
  class(llk) <- "logLik"
  return(llk)
}


