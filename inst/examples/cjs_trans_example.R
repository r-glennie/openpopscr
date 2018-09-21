## Jolly-Seber example 
library(openpopscr)
library(secr)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(lambda0 = 2, sigma = 20, phi = 0.8, sd = 10)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5

# set number of individuals tracked
N <- 100

# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, move = TRUE)



# fit model ---------------------------------------------------------------

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1, 
            sd ~ 1)

start <- list(lambda0 = 2, 
              sigma = 20, 
              phi = 0.5, 
              sd = 2)


obj <- CjsTransientModel$new(par, scrdat, start, num_cores = 4)

obj$par()

obj$calc_llk()

obj$fit()

obj

obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("phi", k = 1)
obj$get_par("sd", k = 1)
