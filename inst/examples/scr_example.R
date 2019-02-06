# SCR example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 2000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh)

# secr fit ----------------------------------------------------------------

screst <- secr.fit(scrdat$capthist(),
                   mask = scrdat$mesh(),
                   detectfn = "HHN")


# openpopscr fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
             sigma  ~ 1)

start <- list(lambda0 = 2, 
              sigma = 20,
              D = 1000)

obj <- ScrModel$new(form, scrdat, start)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("D")
