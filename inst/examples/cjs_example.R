## Jolly-Seber example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(lambda0 = 0.1, sigma = 30, phi = 0.7)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10
primary <- c(rep(1, 3), 
             rep(2, 3), 
             rep(3, 2), 
             rep(4, 2))

# set number of individuals tracked
N <- 1000

# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh)

scrdat2 <- ScrData$new(scrdat$capthist(), mesh, primary = primary)

# fit model ---------------------------------------------------------------

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1)

start <- list(lambda0 = 2, 
              sigma = 20, 
              phi = 0.5)


setThreadOptions(numThreads = 7)

system.time(oo2$calc_pr_capture())

oo <- CjsModel$new(par, scrdat, start)

oo2 <- CjsModel$new(par, scrdat2, start)

oo$par()

oo$calc_llk()

oo$fit()

oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)

