#### SCR density surface example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- secr::make.grid(nx = 7, ny = 7, spacing = 20, detector = "multi")

# make mesh 
mesh <- secr::make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# spatial density effect
x.std <- scale(mesh[,1])
y.std <- scale(mesh[,2])
logD <- log(true_par$D) + 0.2 * x.std - 0.3 * y.std
ihp <- exp(logD) / true_par$D

# set number of occasions to simulate
n_occasions <- 5

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, ihp = ihp)

# plot true density surface
scrdat$plot_mesh(ihp * true_par$D)

# openpopscr fit ----------------------------------------------------------

# formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ x + y)

# starting values 
start <- get_start_values(scrdat)

# create model object 
obj <- ScrModel$new(form, scrdat, start)

# fit model
obj$fit()

# see model results 
obj

# estimated density surface 
Dest <- obj$get_par("D", m = 1:scrdat$n_meshpts())
scrdat$plot_mesh(Dest)

