##### SCR transient model simulation 
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# setup simulations -------------------------------------------------------
set.seed(17417)
nsims <- 100 
ests <- vector(mode = "list", length = nsims)

# set truth 
true_par <- list(D = 1000, lambda0 = 0.5, sigma = 20, sd = 10)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             sd ~ 1, 
             D ~ 1)

# simulation --------------------------------------------------------------
progbar <- utils::txtProgressBar(min = 0, max = nsims, style = 3)  

for (sim in 1:nsims) {
  
  ## progress bar
  Sys.sleep(0.1) 
  utils::setTxtProgressBar(progbar, sim)
  
  # simulate ScrData 
  scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, move = TRUE, print = FALSE)
  
  # get starting values for numerical optimiser  
  start <- list(lambda0 = 0.5, sigma = 20, sd = 10, D = 1000)
  
  # create the model object 
  obj <- ScrTransientModel$new(form, scrdat, start, print = FALSE)
  
  # fit model
  obj$fit()
  
  # store results
  ests[[sim]] <- obj$estimates()$par

}

# extract results 
mu <- exp(sapply(ests, FUN = function(x){x[,1]}))
lcl <- exp(sapply(ests, FUN = function(x){x[,3]}))
ucl <- exp(sapply(ests, FUN = function(x){x[,4]}))
# distributions 
summary(t(mu))
# confidence interval coverage 
sum(lcl[1,] < true_par$lambda0 & ucl[1,] > true_par$lambda0)
sum(lcl[2,] < true_par$sigma & ucl[2,] > true_par$sigma)
sum(lcl[3,] < true_par$sd & ucl[3,] > true_par$sd)
sum(lcl[4,] < true_par$D & ucl[4,] > true_par$D)
