# autoini function

#' Get initial values
#'
#' @param obj ScrData
#' @param model model class
#'
#' @return initial values for parameters 
#' @export
get_start_values <- function(obj, model = "ScrModel") {
  auto.sigma <- RPSV(obj$capthist(), CC = TRUE)
  est_encrate <- obj$encrate() 
  # compute unit encounter rate for given sigma
  r <- obj$distances()
  encrate <- colSums(exp(-r^2 / (2 * auto.sigma^2)))
  unit_encrate <- mean(encrate)
  # guess lambda0 
  auto <- list(lambda0 = est_encrate / unit_encrate) 
  auto$sigma <- auto.sigma
  # compute detection probability given lambda0, sigma
  pdot <- mean(1 - exp(-auto$lambda0 * encrate * obj$n_occasions()))
  auto$D <- obj$n() / (obj$area() * pdot)
  if (model %in% c("ScrTransientModel", "CjsTransientModel", "JsTransientModel")) {
    auto$sd <- obj$encrange()^2 - mean(obj$encrange(each = TRUE))^2
    if (auto$sd < 0) {
      auto$sd <- 1e-10
    } else {
      auto$sd <- sqrt(auto$sd)
    }
  } 
  if (model %in% c("CjsModel", "CjsTransientModel", "JsModel", "JsTransientModel")) {
    range_seen <- sapply(1:obj$n(), FUN = function(i) {
      redcap <- rowSums(obj$capthist()[i,,])
      seen <- which(redcap > 0)
      first <- min(seen)
      last <- max(seen)
      return(c(first, last)) 
    })
    diff <- diff(range_seen)
    diff <- diff[diff > 0]
    auto$phi <- 1 - 1 / mean(diff)
    auto$phi <- auto$phi ^ (1 / mean(diff(obj$time())))
  }
  if (model %in% c("JsModel", "JsTransientModel")) {
    auto$beta <- sum(range_seen[1,]==1)/ncol(range_seen)
  }
  if (model %in% c("CjsModel", "CjsTransientModel")) auto$D <- NULL
  if (any(sapply(auto, is.na))) stop("autoini produces NA, set starting values manually.")
  return(auto)
}

