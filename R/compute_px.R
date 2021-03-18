# estimate posterior probability of activity centre on mesh 
# for non-moving activity centres 
#' Title
#'
#' @param mod fitted non-transient model 
#'
#' @return matrix of posterior probabilities mesh point by individual 
#' @export
compute_px <- function(mod) {
  # check model type
  if (!(class(mod)[1] %in% c("JsModel", "ScrModel"))) stop("compute_px only works for JsModel and ScrModel types")
  # initial distribution 
  pr0 <- mod$calc_initial_distribution()
  # compute probability of capture histories 
  # across all individuals, occasions and traps 
  pr_capture <- mod$calc_pr_capture()
  # compute lalpha for each individual
  n <- mod$data()$n()
  n_occasions <- mod$data()$n_occasions()
  n_meshpts <- mod$data()$n_meshpts() 
  # get tpms for state model 
  nstates <- mod$nstates()
  tpms <- mod$calc_tpms()
  
  lalpha <- openpopscr:::C_calc_alpha(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, mod$data()$n()))
  lbeta <- openpopscr:::C_calc_beta(n, n_occasions, n_meshpts, pr0, pr_capture, tpms, nstates, rep(0, mod$data()$n()))
  px <- matrix(0, nr = mod$data()$n_meshpts(), nc = mod$data()$n())
  for (i in 1:length(lalpha)) {
    pi <- lalpha[[i]][,,mod$data()$n_occasions(), drop=FALSE]
    pmax <- max(pi)
    psum <- log(sum(exp(pi - pmax))) + pmax
    pxs <- exp(pi - psum) 
    px[,i] <- rowSums(pxs)
  }
  return(px)
}
