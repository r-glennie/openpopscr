# SCR example 
library(openpopscr)
library(secr)


# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# detector factor covariate
type <- sample(0:1, size = nrow(detectors), replace = TRUE)
covariates(detectors) <- data.frame(type = type)

# detector continuous covariate 
cover <- runif(nrow(detectors))
covariates(detectors)$cover <- cover

# temporal factor covariate 
temp <- factor(c(0, 0, 1, 1, 0))

# temporal continuous covariate
temp2 <- runif(5) * 10

# simulate ScrData 
scrdat <- simulate_scr(true_par, 
                       n_occasions, 
                       detectors, 
                       mesh, 
                       seed = 15483)

# secr fit ----------------------------------------------------------------

screst <- secr.fit(scrdat$capthist(),
                   mask = scrdat$mesh(),
                   model = list(sigma ~ cover + temp, lambda0 ~ type + temp2), 
                   timecov = data.frame(temp = temp, temp2 = temp2), 
                   detectfn = "HHN")


# openpopscr fit ----------------------------------------------------------

scrdat <- ScrData$new(scrdat$capthist(), 
                      scrdat$mesh(), 
                      scrdat$time(), 
                      list(type = type, cover = cover, temp = temp, temp2 = temp2),
                      c("j", "j", "k", "k"))

form <- list(lambda0 ~ type + temp2, 
             sigma  ~ cover + temp)

start <- list(lambda0 = 2, 
              sigma = 20,
              D = 1000)

obj <- ScrModel$new(form, scrdat, start, num_cores = 4)

# compute initial likelihood 
#obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1, j = 1)
obj$get_par("sigma", k = 1, j = 1)
obj$get_par("D")
