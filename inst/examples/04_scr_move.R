# SCR transient model example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 7)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 0.5, sigma = 20, sd = 10)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, move = TRUE, seed = 15483)

# stationary fit ----------------------------------------------------------------

stat <- ScrModel$new(list(lambda0 ~ 1, sigma ~ 1, D ~ 1), 
                     scrdat, 
                     list(lambda0 = 0.5, sigma = 20, D = 1000))
stat$fit()
stat

# transient fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             sd ~ 1, 
             D ~ 1)

# get automatic start values 
start <- get_start_values(scrdat, model = "ScrTransientModel")
# update start values using stationary model fit 
start$lambda0 <- exp(stat$par()$lambda0[[1]])
start$sigma <- exp(stat$par()$sigma[[1]])
start$D <- exp(stat$par()$D[[1]])

# create model object
obj <- ScrTransientModel$new(form, scrdat, start)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# compare stationary and transience models 
AIC(stat, obj)
