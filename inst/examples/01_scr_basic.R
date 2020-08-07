#### Basic SCR example
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------
# set truth 
true_par <- list(D = 1000, lambda0 = 0.2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "multi")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh)

# openpopscr fit ----------------------------------------------------------

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             D ~ 1)

# get starting values for numerical optimiser  
start <- list(lambda0 = 0.2, sigma = 20, D = 1000)

# create the model object 
obj <- ScrModel$new(form, scrdat, start)
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1, j = 1)
obj$get_par("sigma", k = 1, j = 1)
obj$get_par("D")
