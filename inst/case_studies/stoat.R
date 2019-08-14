## Analysis of Stoat DNA data 
## Data avaliable from secr R package 
## Reference: Gleeson, D. M., Byrom, A. E. and Howitt, R. L. J. (2010) 
##            Non-invasive methods for genotyping of stoats (Mustela erminea) in New Zealand: potential for field applications. 
##            New Zealand Journal of Ecology 34, 356–359. Available on-line at http://www.newzealandecology.org.

library(openpopscr)
library(secr)

# data --------------------------------------------------------------------

# rough activity range estimate
emp_sigma <- RPSV(stoatCH)

# create mesh with buffer ~3*emp_sigma 
mesh <- make.mask(traps(stoatCH), buffer = 1000, type = "trapbuffer")

# create data object
dat <- ScrData$new(stoatCH, mesh)

# output data 
dat


# fit basic model ---------------------------------------------------------

# get starting values
start <- get_start_values(dat)

# setup formulae
f0 <- list(lambda0 ~ 1, sigma ~ 1, D ~ 1)

# create basic model
m0 <- ScrModel$new(f0, dat, start)

# fit model
m0$fit()

# update starting values 
ests <- m0$estimates()$par[,1]
start2 <- list(lambda0 = exp(ests[1]), 
              sigma = exp(ests[2]), 
              D = exp(ests[3]))

# detection models --------------------------------------------------------

# add time 
flamt <- list(lambda0 ~ s(t, k = 3), sigma ~ 1, D ~ 1)
mlamt <- ScrModel$new(flamt, dat, start2)
mlamt$fit()

fsigt <- list(lambda0 ~ 1, sigma ~ s(t, k = 3), D ~ 1)
msigt <- ScrModel$new(fsigt, dat, start2)
msigt$fit()

AIC(m0, mlamt, msigt)
# no evidence of a time effect 

# spatial models ----------------------------------------------------------

fDx <- list(lambda0 ~ 1, sigma ~ 1, D ~ s(x, k = 3))
mDx <- ScrModel$new(fDx, dat, start2)
mDx$fit()

fDy <- list(lambda0 ~ 1, sigma ~ 1, D ~ s(y, k = 3))
mDy <- ScrModel$new(fDy, dat, start2)
mDy$fit()

fDxy <- list(lambda0 ~ 1, sigma ~ 1, D ~ s(x, y, k = 5))
mDxy <- ScrModel$new(fDxy, dat, start2)
mDxy$fit()

fDx5 <- list(lambda0 ~ 1, sigma ~ 1, D ~ s(x, k = 5))
mDx5 <- ScrModel$new(fDx5, dat, start2)
mDx5$fit()

AIC(m0, mDx, mDy, mDxy, mDx5)
# Very weak support for mDx, decided to stick with m0


# inference ---------------------------------------------------------------

# look at estimates 
m0 
m0$get_par("D")
m0$predict("D", se = TRUE)[1,]

# plot detection fn 
est <- exp(m0$estimates()$par[1:2, 1])
m0$detectfn()$plot(est, xlim = c(0, 500), ylim = c(0, 0.06))

# plot hazard fn 
m0$detectfn()$plot(est, h = TRUE, xlim = c(0, 500))


