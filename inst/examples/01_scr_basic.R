#### Basic SCR example
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 3)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh)

# openpopscr fit ----------------------------------------------------------

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1)

# get starting values for numerical optimiser  
start <- get_start_values(scrdat)

# create the model object 
obj <- ScrModel$new(form, scrdat, start)

# compute initial likelihood to see if start is reasonable
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1, j = 1)
obj$get_par("sigma", k = 1, j = 1)
obj$get_par("D")
