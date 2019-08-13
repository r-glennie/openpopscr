## Jolly-Seber example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20, phi = 0.5, beta = 0.25, sd = 10)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5

# simulate ScrData 
scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, move = TRUE, seed = 2952)

# fit model ---------------------------------------------------------------

# create formulae 
par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            beta ~ 1, 
            phi ~ 1, 
            sd ~ 1, 
            D ~ 1)

# get starting values 
start <- get_start_values(scrdat, model = "JsTransientModel")

# create model object 
oo <- JsTransientModel$new(par, scrdat, start)

# compute initial likelihood
oo$calc_llk()

# fit model 
oo$fit()

# see results 
oo

oo$get_par("lambda0", k = 1, j = 1)
oo$get_par("sigma", k = 1, j = 1)
oo$get_par("phi", k = 1)
oo$get_par("beta", k = 1)
oo$get_par("sd", k = 1)
oo$get_par("D")
