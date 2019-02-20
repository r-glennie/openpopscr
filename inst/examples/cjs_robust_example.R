## CJS Robust Design example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(lambda0 = 2, sigma = 30, phi = 0.7)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "multi")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of total occasions to simulate
n_occasions <- 10

# set primary periods 
primary <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2))

# set number of individuals tracked
N <- 800

# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, primary = primary)


# fit model ---------------------------------------------------------------

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ primary)

start <- get_start_values(scrdat, model = "CjsModel")

oo <- CjsModel$new(par, scrdat, start)

oo$par()

oo$calc_llk()

oo$fit()

oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)


# test openCR 

spltcap <- split(scrdat$capthist(), as.factor(primary), byoccasion = TRUE)

fit.nsp <- openCR.fit(spltcap, type = "CJS")

fit.sp <- openCR.fit(spltcap, type = "CJSsecr", mask = scrdat$mesh(), trace = TRUE)
fit.sp
