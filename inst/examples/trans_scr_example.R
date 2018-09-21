# SCR example 
library(openpopscr)
library(secr)


# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20, sd = 10)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, move = TRUE, seed = 15483)

# stationary fit ----------------------------------------------------------------

stat <- ScrModel$new(list(lambda0 ~ 1, sigma ~ 1), 
                     scrdat, 
                     list(lambda0 = 2, sigma = 20, D = 1000), 
                     num_cores = 4)

stat$fit()

stat$get_par("lambda0", k = 1)
stat$get_par("sigma", k = 1)
stat$get_par("D")

# transient fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             sd ~ 1)

start <- list(lambda0 = 2, 
              sigma = 20,
              sd = 10, 
              D = 1000)

obj <- ScrTransientModel$new(form, scrdat, start, num_cores = 4)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("sd", k = 1)
obj$get_par("D")



