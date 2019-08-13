## JS Robust Design example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 3)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 0.1, sigma = 40, phi = 0.7, beta = 0.6)

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

# create formulae 
par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1, 
            beta ~ 1,
            D ~ 1)

# get start values
start <- get_start_values(scrdat, model = "JsModel")

# create model object 
oo <- JsModel$new(par, scrdat, start)

# compute initial likelihood 
oo$calc_llk()

# fit model 
oo$fit()

# see results 
oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)
oo$get_par("beta", k = 1)
oo$get_par("D")
