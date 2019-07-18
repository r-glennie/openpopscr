### Stratfied analysis example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 3)

# simulate data
# set true parameters 
true_par <- list(s1 = list(D = 1000, lambda0 = 2, sigma = 20), 
                 s2= list(D = 500, lambda0 = 2, sigma = 20), 
                 s3 = list(D = 1000, lambda0 = 2, sigma = 10))

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "multi")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate for each stratum
scrdat1 <- simulate_scr(true_par[[1]], n_occasions, detectors, mesh, seed = 15483)
scrdat2 <- simulate_scr(true_par[[2]], n_occasions, detectors, mesh, seed = 79696)
scrdat3 <- simulate_scr(true_par[[3]], n_occasions, detectors, mesh, seed = 43523)

scrdat <- list(s1 = scrdat1, s2 = scrdat2, s3 = scrdat3)


# create formula for shared parameters 
shared_form <- list(lambda0 ~ 1)

# create formulae for private parameters 
private_form <- list(s1 = list(sigma ~ 1),
                     s2 = list(sigma ~ 1), 
                     s3 = list(sigma ~ 1))


# set starting values 
start <- list(lambda0 = 2, sigma = 20, D = 1000)

# create model object 
obj <- StrataModel$new(scrdat, "ScrModel", shared_form, private_form, start)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see results 
obj 


