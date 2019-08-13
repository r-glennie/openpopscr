#### SCR with covariates on detection example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------
## simulation functions in openpopscr cannot simulate with covariate effects
## here, we simulate with constant detection parameters, but fit models 
## with covariates, just to see how it would be done for real data. 

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrdat <- simulate_scr(true_par, 
                       n_occasions, 
                       detectors, 
                       mesh, 
                       seed = 15483)

# add detector factor covariate
type <- factor(sample(0:1, size = nrow(detectors), replace = TRUE))
scrdat$add_covariate("type", type, "j")

# detector continuous covariate 
cover <- runif(nrow(detectors))
scrdat$add_covariate("cover", cover, "j")

# temporal factor covariate 
temp <- factor(c(0, 0, 1, 1, 0))
scrdat$add_covariate("temp", temp, "k")

# temporal continuous covariate
temp2 <- runif(5) * 10
scrdat$add_covariate("temp2", temp2, "k")


# openpopscr fit ----------------------------------------------------------
# have encounter rate depend on detector type and temporal covariate 2 
# have encounter range depend on vegetative cover around detector and temporal covariate 1
form <- list(lambda0 ~ type + temp2, 
             sigma  ~ cover + temp, 
             D ~ 1)

start <- get_start_values(scrdat)

obj <- ScrModel$new(form, scrdat, start)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj
