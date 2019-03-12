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
################################################################################
# openpopscr project: open population spatial capture-recapture 
#
# file description: encounter function / detection function object 
#
################################################################################

#' Detection function 
#' 
#' @description Object represents the encounter function / detection function. 
#' \itemize{
#'   \item parameters: a character vector of parameter names
#'   \item fn: either a character denoting a encoutner rate function / detection function or a function
#'   that takes parameter vector and distance and returns encounter rate or probability 
#'   \item prob: if fn supplied by user return probability then set this to TRUE, otherwise it is FALSE 
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item g: returns detection probability at a given distance and parameter values 
#'  \item h: returns encounter rate at a given distance and parameter values 
#'  \item plot: plots detection function if prob = TRUE, otherwise encounter function 
#'  \item print: prints out basic encounter information 
#' }
#' 
DetFn <- R6Class("DetFn", 
  public = list( 
    initialize = function(parameters = NULL, 
                          fn = NULL, 
                          prob = FALSE, 
                          link2response = NULL, 
                          response2link = NULL) {
      if (is.null(fn)) fn <- "HHN"
      if (class(fn) == "function") {
        private$fn_ <- fn
        private$parameters_ <- parameters 
        private$prob_ <- prob 
        private$link2response_ <- link2response
        private$response2link_ <- response2link 
      } else if (class(fn) == "character") {
        if (fn == "HN") {
          private$fn_ <- function(x, par) {par[[1]]*exp(-(x^2)/(2*par[[2]]^2))}
          private$prob_ <- TRUE
          private$parameters_ <- c("g0", "sigma")
          private$link2response_ <- list("plogis", "exp")
          private$response2link_ <- list("qlogis", "log")
        } else if (fn == "HHN") {
          private$fn_ <- function(x, par) {par[[1]]*exp(-(x^2)/(2*par[[2]]^2))}
          private$prob_ <- FALSE
          private$parameters_ <- c("lambda0", "sigma")
          private$link2response_ <- list("exp", "exp")
          private$response2link_ <- list("log", "log")
        }
      } 
    },
    
    print = function(i = ".") {
       print(private$fn_)
    }, 
    
    g = function(x, par) {
      if (private$prob_) return(private$fn_(x, par))
      return(1 - exp(-private$fn_(x, par)))
    }, 
    
    h = function(x, par) {
      if (!private$prob_) return(private$fn_(x, par))
      return(-log(1 - private$fn_(x, par)))
    },
    
    pars = function() {
      return(private$parameters_)
    }, 
    
    par = function(i) {
      return(private$parameters_[i])
    }, 
    
    npars = function() {
      return(length(private$parameters_))
    }, 
    
    link2response = function(i) {
      return(private$link2response_[i])
    },
    
    response2link = function(i) {
      return(private$response2link_[i])
    }, 
    
    plot = function(par, 
                    h = FALSE, 
                    xlim = c(0,1), 
                    ylim = NULL, 
                    xname = "Distance", 
                    yname = NULL,
                    main = NULL,
                    ...) {
      gfn <- function(x) {self$g(x, par)}
      hfn <- function(x) {self$h(x, par)}
      if (!h) {
         if (is.null(yname)) yname <- "Probability"
         if (is.null(main)) main <- "Detection Function"
         if (is.null(ylim)) ylim <- c(0, 1)
         curve(gfn, from = xlim[1], to = xlim[2], ylim = ylim, xname = xname, ylab = yname, main = main, ...)
      }
      if (h) {
        if (is.null(yname)) yname <- "Encounter Rate"
        if (is.null(main)) main <- "Encounter Function"
        curve(hfn, from = xlim[1], to = xlim[2], ylim = ylim, xname = xname, ylab = yname, main = main, ...)
      }
      invisible(list(par, h, xlim, ylim, xname, yname, main, ...))
    }
    
  ), 
  
  private = list(
    parameters_ = NULL, 
    fn_ = NULL, 
    prob_ = NULL, 
    link2response_ = NULL, 
    response2link_ = NULL
  )
)








