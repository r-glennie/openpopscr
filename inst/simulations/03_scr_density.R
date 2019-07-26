#### SCR density surface simulation 
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# setup simulations -------------------------------------------------------

set.seed(14919)
nsims <- 100 
ests <- vector(mode = "list", length = nsims)

# set truth 
true_par <- list(D = 1000, lambda0 = 0.2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# spatial density effect
x.std <- scale(mesh[,1])
y.std <- scale(mesh[,2])
logD <- log(true_par$D) + 0.2 * x.std - 0.3 * y.std
ihp <- exp(logD) / true_par$D

# set number of occasions to simulate
n_occasions <- 5 

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ x + y)

# simulation --------------------------------------------------------------
progbar <- utils::txtProgressBar(min = 0, max = nsims, style = 3)  

for (sim in 1:nsims) {
  
  ## progress bar
  Sys.sleep(0.1) 
  utils::setTxtProgressBar(progbar, sim)
  
  # simulate ScrData 
  scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, ihp = ihp, print = FALSE)
  
  # get starting values for numerical optimiser  
  start <- get_start_values(scrdat)
  
  # create the model object 
  obj <- ScrModel$new(form, scrdat, start, print = FALSE)
  
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
sum(lcl[1,] < log(true_par$lambda0) & ucl[1,] > log(true_par$lambda0), na.rm=TRUE)
sum(lcl[2,] < log(true_par$sigma) & ucl[2,] >log(true_par$sigma),na.rm=TRUE)
sum(lcl[3,] < log(true_par$D) & ucl[3,] > log(true_par$D), na.rm=TRUE)
sum(lcl[4,] < 0.2 & ucl[4,] > 0.2, na.rm=TRUE)
sum(lcl[5,] < -0.3 & ucl[5,] > -0.3, na.rm=TRUE)

