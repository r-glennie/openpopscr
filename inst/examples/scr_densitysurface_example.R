# SCR density surface example 
library(openpopscr)
library(secr)
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
logD <- log(true_par$D) + 0.2 * x.std
ihp <- exp(logD) / true_par$D

# set number of occasions to simulate
n_occasions <- 5

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, ihp = ihp)

# secr fit ----------------------------------------------------------------

screst <- secr.fit(scrdat$capthist(),
                   mask = scrdat$mesh(),
                   model = list(D ~ x), 
                   detectfn = "HHN")


# openpopscr fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ 1)

start <- get_start_values(scrdat, model = "ScrTransientModel")

obj <- ScrModel$new(form, scrdat, start)

# fit model
obj$fit()

fit <- openCR.fit(scrdat$capthist(), mask = scrdat$mesh(), type = "secrD", trace = TRUE)

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("sd", k = 1)
obj$get_par("D")
