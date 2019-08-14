#### Jolly-Seber SCR density surface example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 200, lambda0 = 0.5, sigma = 20, phi = 0.8, beta = 0.5)

# make detectors array 
detectors <- secr::make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- secr::make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# spatial density effect
x.std <- scale(mesh[,1])
y.std <- scale(mesh[,2])
logD <- log(true_par$D) +  2 * x.std - 3 * y.std
ihp <- exp(logD) / true_par$D

# set number of occasions to simulate
n_occasions <- 5

# simulate ScrData 
scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, ihp = ihp, seed = 1481)

# openpopscr fit ----------------------------------------------------------

# create formulae 
form <- list(lambda0 ~ 1, 
            sigma ~ 1,
            beta ~ 1, 
            phi ~ 1,
            D ~ x + y)

# get starting values 
start <- get_start_values(scrdat, model = "JsModel")

# create model object 
obj <- JsModel$new(form, scrdat, start)

# compute initial likelihood
obj$calc_llk()

# fit model
obj$fit()

# see results 
obj

# plot density surface
D <- obj$get_par("D", m = 1:scrdat$n_meshpts()) 
scrdat$plot_mesh(D)
