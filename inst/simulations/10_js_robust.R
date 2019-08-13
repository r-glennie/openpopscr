#### Robust JS simulation
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 30)

# setup simulations -------------------------------------------------------

set.seed(184821)
nsims <- 100 
ests <- vector(mode = "list", length = nsims)

# set truth 
true_par <- list(D = 1000, lambda0 = 0.5, sigma = 20, phi = 0.8, beta = 0.5)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of total occasions to simulate
n_occasions <- 10

# set primary periods 
primary <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2))

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma ~ 1,
             beta ~ 1, 
             phi ~ 1,
             D ~ x + y)

# simulation --------------------------------------------------------------
progbar <- utils::txtProgressBar(min = 0, max = nsims, style = 3)  

for (sim in 1:nsims) {
  
  ## progress bar
  Sys.sleep(0.1) 
  utils::setTxtProgressBar(progbar, sim)
  
  # simulate ScrData 
  scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, primary = primary, print = FALSE)
  
  # get starting values for numerical optimiser  
  start <- get_start_values(scrdat, model = "JsModel")
  
  # create the model object 
  obj <- JsModel$new(form, scrdat, start, print = FALSE)
  
  # fit model
  obj$fit()
  
  # store results
  ests[[sim]] <- obj$estimates()$par

}

# extract results 
mu <- sapply(ests, FUN = function(x){x[,1]})
lcl <- sapply(ests, FUN = function(x){x[,3]})
ucl <- sapply(ests, FUN = function(x){x[,4]})
# distributions 
summary(t(mu))
# confidence interval coverage 
sum(lcl[1,] < log(true_par$lambda0) & ucl[1,] > log(true_par$lambda0))
sum(lcl[2,] < log(true_par$sigma) & ucl[2,] > log(true_par$sigma))
sum(lcl[3,] < qlogis(true_par$phi) & ucl[3,] > qlogis(true_par$phi))
b <- (1 - true_par$beta) / (n_occasions - 1)
b <- log(b/true_par$beta)
sum(lcl[4,] < b & ucl[4,] > b)
sum(lcl[5,] < log(true_par$D) & ucl[5,] > log(true_par$D))
