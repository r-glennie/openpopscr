## JS Robust Design example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 1, sigma = 40, phi = 0.7, beta = 0.6)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 200, nx = 64, ny = 64, type = "trapbuffer")

# set number of total occasions to simulate
n_occasions <- 10

# set primary periods 
primary <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2))

# simulate ScrData 
scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, primary = primary)


# fit model ---------------------------------------------------------------

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1, 
            beta ~ primary)

start <- get_start_values(scrdat, model = "JsModel")


oo <- JsModel$new(par, scrdat, start)

oo$par()

oo$calc_llk()

oo$fit()

oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)
oo$get_par("beta", k = 1)
oo$get_par("D")
