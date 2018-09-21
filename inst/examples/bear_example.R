# bear example 
library(openpopscr)

beardata <- readRDS("bear_data.Rds")
capthist <- beardata$capthist 
mesh <- beardata$mesh

bear_dat <- ScrData$new(capthist, mesh)

RPSV(capthist)
autoini(capthist, mesh)

form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             sd ~ 1)

start <- list(lambda0 = 0.2, 
              sigma = 2000,
              sd = 1000, 
              D = 0.15)

obj <- ScrTransientModel$new(form, bear_dat, start, num_cores = 4)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1)
obj$get_par("sigma", k = 1)
obj$get_par("sd", k = 1)
obj$get_par("D")


