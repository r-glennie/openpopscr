# SCR transect example 
library(openpopscr)
library(secr)


# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 100, lambda0 = 0.2, sigma = 400)

# make transects 
starts <- seq(100, 1000, 200)
transect_list <- vector(mode = "list", length = 5)
for (i in 1:5) {
  transect_list[[i]] <- data.frame(x = starts[i], y = seq(0, 1000, 200))
}
detectors <- make.transect(transect_list)

# make mesh 
mesh <- make.mask(detectors, buffer = 500)

# set number of occasions to simulate
n_occasions <- 10 

# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, seed = 15483)

# secr fit ----------------------------------------------------------------

screst <- secr.fit(scrdat$capthist(),
                   mask = scrdat$mesh(),
                   detectfn = "HHN")


# openpopscr fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
             sigma  ~ 1)

start <- list(lambda0 = 1, 
              sigma = 20,
              D = 100)

obj <- ScrModel$new(form, scrdat, start, num_cores = 4)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("D")
