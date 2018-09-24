## Data 

#' Stoat data as included in secr package 
#'
#' This is a copy of the data included in the secr package. 
#' Data of A. E. Byrom from a study of stoats (Mustela erminea) in New Zealand. Individuals were identified from DNA in hair samples.
#'
#' @format A list of two data frames 
#' \describe{
#'   \item{captures}{Session, ID, Occasion, and Detector} 
#'   \item{traps}{Detector, x, and y}
#' }
"stoat"

#' Fitted models from ScrModel vignette 
#'
#' Includes all fitted models from vignette. 
#'
#' @format A list of two data frames 
#' \describe{
#'   \item{mod}{standard SCR model} 
#'   \item{mod_detage}{SCR model with detector covariate}
#'   \item{stat}{stationary SCR model to transient data}
#'   \item{trans}{transient SCR model to transient data}
#' }
"scrmodel_example"

#' Fitted models from CjsModel vignette 
#'
#' Includes all fitted models from vignette. 
#'
#' @format A list of two data frames 
#' \describe{
#'   \item{mod}{standard SCR model} 
#'   \item{mod_surv}{SCR model with detector covariate}
#'   \item{stat}{stationary SCR model to transient data}
#'   \item{trans}{transient SCR model to transient data}
#' }
"cjsmodel_example"

#' Fitted models from JsModel vignette 
#'
#' Includes all fitted models from vignette. 
#'
#' @format A list of two data frames 
#' \describe{
#'   \item{mod}{standard SCR model} 
#'   \item{mod_seas}{SCR model with detector covariate}
#'   \item{stat}{stationary SCR model to transient data}
#'   \item{trans}{transient SCR model to transient data}
#' }
"jsmodel_example"


#' Fitted models from Strata model vignette 
#'
#' Includes all fitted models from vignette. 
#'
#' @format A fitted model object  
#' \describe{
#'   \item{obj}{Strata model object}
#' }
"stratamodel_example"
