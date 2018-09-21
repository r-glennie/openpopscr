library(openpopscr)
library(secr)


# SCR ---------------------------------------------------------------------

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


# create model object

shared_form <- list(lambda0 ~ 1)

private_form <- list(s1 = list(sigma ~ 1),
                     s2 = list(sigma ~ 1), 
                     s3 = list(sigma ~ 1))


start <- list(lambda0 = 2, sigma = 20, D = 1000)

obj <- StrataModel$new(scrdat, "ScrModel", shared_form, private_form, start)



# fit model

inillk <- obj$calc_llk()

obj$fit()


# CJS ---------------------------------------------------------------------

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(s1 = list(lambda0 = 2, sigma = 20, phi = 0.5),
                 s2 = list(lambda0 = 2, sigma = 20, phi = 0.7), 
                 s3 = list(lambda0 = 2, sigma = 20, phi = 0.9))

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5

# set number of individuals tracked
N <- 100

# simulate ScrData 
scrdat1 <- simulate_cjs_openscr(true_par[[1]], N, n_occasions, detectors, mesh, seed = 12951)
scrdat2 <- simulate_cjs_openscr(true_par[[2]], N, n_occasions, detectors, mesh, seed = 25825)
scrdat3 <- simulate_cjs_openscr(true_par[[3]], N, n_occasions, detectors, mesh, seed = 51811)

scrdat <- list(s1 = scrdat1, s2 = scrdat2, s3 = scrdat3)

# fit model ---------------------------------------------------------------

shared_form <- list(lambda0 ~ 1, sigma ~ 1)

private_form <- list(s1 = list(phi ~ 1),
                     s2 = list(phi ~ 1), 
                     s3 = list(phi ~ 1))


start <- list(lambda0 = 2, sigma = 20, phi = 0.6, D = 1000)

obj <- StrataModel$new(scrdat, "CjsModel", shared_form, private_form, start)

obj$calc_llk()

obj$fit()

obj$get_object(3)$get_par("phi", k =1)

# CJS with transience ----------------------------------------------------------

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(s1 = list(lambda0 = 2, sigma = 20, phi = 0.5, sd = 10),
                 s2 = list(lambda0 = 2, sigma = 20, phi = 0.7, sd = 15), 
                 s3 = list(lambda0 = 2, sigma = 20, phi = 0.9, sd = 30))

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5

# set number of individuals tracked
N <- 100

# simulate ScrData 
scrdat1 <- simulate_cjs_openscr(true_par[[1]], N, n_occasions, detectors, mesh, move = TRUE, seed = 12941)
scrdat2 <- simulate_cjs_openscr(true_par[[2]], N, n_occasions, detectors, mesh, move = TRUE, seed = 25835)
scrdat3 <- simulate_cjs_openscr(true_par[[3]], N, n_occasions, detectors, mesh, move = TRUE, seed = 51821)

scrdat <- list(s1 = scrdat1, s2 = scrdat2, s3 = scrdat3)

# fit model ---------------------------------------------------------------

shared_form <- list(lambda0 ~ 1, sigma ~ 1)

private_form <- list(s1 = list(phi ~ 1, sd ~ 1),
                     s2 = list(phi ~ 1, sd ~ 1), 
                     s3 = list(phi ~ 1, sd ~ 1))


start <- list(lambda0 = 2, sigma = 20, phi = 0.6, sd = 10, D = 1000)

obj <- StrataModel$new(scrdat, "CjsTransientModel", shared_form, private_form, start)

obj$calc_llk()

obj$fit()

obj$get_object(3)$get_par("phi", k =1)
obj$get_object(3)$get_par("sd", k =1)
