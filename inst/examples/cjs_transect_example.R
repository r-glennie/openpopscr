## Jolly-Seber example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(lambda0 = 1, sigma = 400, phi = 0.7)

# make transects 
starts <- seq(100, 1000, 200)
transect_list <- vector(mode = "list", length = 5)
for (i in 1:5) {
  transect_list[[i]] <- data.frame(x = starts[i], y = seq(0, 1000, 200))
}
detectors <- make.transect(transect_list)

# make mesh 
mesh <- make.mask(detectors, buffer = 500, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10

# set number of individuals tracked
N <- 200

# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, seed = 12952)

# fit model ---------------------------------------------------------------

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1)

start <- get_start_values(scrdat, model = "CjsModel")


oo <- CjsModel$new(par, scrdat, start)

oo$par()

oo$calc_llk()

oo$fit()

oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)
