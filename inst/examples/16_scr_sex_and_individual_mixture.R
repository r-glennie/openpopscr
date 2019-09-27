#### Basic partially observed mixture SCR example
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 6)

# simulate data -----------------------------------------------------------
set.seed(53919)

## Set parameters 
D <- 5000
lambda0 <- rep(0.1, 4) # one for each state 
sigma <-  exp(log(20) + c(0, 1, -0.5, 0.5)) # one for each state 

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
K <- 10

## Simulate activity centres
A <- nrow(mesh) * attr(mesh, "area") / 100
N <- rpois(1, D*A)
pt <- sample(1:nrow(mesh), size = N, replace = TRUE)
x <- mesh[pt, 1]
y <- mesh[pt, 2]

## Simulate mixture 
nstates <- 4
delta <- c(c(0.3, 0.7) * 0.6, c(0.4, 0.6) * 0.4) # state membership probability 
mix <- sample(1:4, size = N, replace = TRUE, prob = delta)
obsmix <- rep(NA, length(mix))
cmix <- NULL

## Simulate survey 
cap <- data.frame(session = numeric(), 
                  ID = numeric(), 
                  occasion = numeric(), 
                  trap = numeric())

seen <- rep(FALSE, N)
id <- rep(0, N)
for (k in 1:K) {
  for (i in 1:N) {
    d2 <- (x[i] - detectors[,1])^2 + (y[i] - detectors[,2])^2
    er <- lambda0[mix[i]] * exp(-d2 / (2 * sigma[mix[i]]^2))
    c <- rpois(length(er), er)
    if (any(c > 0)) {
      if (!seen[i]) {
        id[i] <- max(id) + 1
        if (is.na(obsmix[i])) cmix <- c(cmix, mix[i])
        obsmix[i] <- mix[i]
        seen[i] <- TRUE
      } 
      dets <- which(c > 0)
      for (r in 1:length(dets)) {
        nc <- c[dets[r]]
        rec <- data.frame(session = rep(1, nc), 
                          ID = rep(id[i], nc), 
                          occasion = rep(k, nc), 
                          trap = rep(dets[r], nc))
        cap <- rbind(cap, rec)
      }
    }
  }
}
if (max(cap$occasion) != K) cap <- rbind(cap, data.frame(session = 1, ID = "NONE", occasion = K, trap = 1))
ch <- make.capthist(cap, detectors)

# create data object 
scrdat <- ScrData$new(ch, mesh = mesh)

# create observed sex covariate
si <- rep(NA, scrdat$n()) # NA for those individuals where sex is unknown 
si[1:250] <- cmix[1:250]
si[si %in% c(1, 2)] <- 1 
si[si %in% c(3, 4)] <- 2 
# add to data 
scrdat$add_covariate("sex", si, "i")

# openpopscr fit ----------------------------------------------------------

## create state model 
# Here, we want basic mixture model, no state-switching so the structure is 
# four independent states. "0" in struct indicates no switches. "." is always
# on the diagonal. 
struct <- matrix("0", nr = 4, nc = 4)
diag(struct) <- "."

# Here, I create four states male1, male2, female1, female2 so that
# each sex has two latent sub-populations. 
# I must specify an initial state membership (I just made it uniform) and a tpm (no switching)
# cov is a dataframe that maps the states to any corresponding individual level covariates, so here
# I specify that sex is 1 (male) for states male1, male2 and is 2 for the female states. I also
# specify a subpop covariate. 
statemod <- StateModel$new(data = scrdat, 
                           names = c("male1", "male2", "female1", "female2"), 
                           structure = struct, 
                           start = list(delta = c(0.25, 0.25, 0.25, 0.25), tpm = diag(4)), 
                           cov = data.frame(sex = factor(c(1, 1, 2, 2)), 
                                            subpop = factor(c(1, 2, 1, 2))))

## create formulae 
# For sigma, I want a sex-specific effect and a subpop-specific effect given sex. 
# If I had not created sex, subpop in cov above, then there is only one built-in variable
# I can use for mixtures: state. "state" is equivalent to a factor variable, one level for each
# hidden state. There is no checking whether the model you are trying to fit makes sense. 
# For example sigma ~ sex*subpop + state is nonsense because sex*subpop is equivalent to state and 
# so your model has redundant parameters. 
form <- list(lambda0 ~ 1, 
             sigma  ~ sex + subpop,
             D ~ 1)

# get starting values for numerical optimiser  
start <- get_start_values(scrdat)

# create the model object 
obj <- ScrModel$new(form, scrdat, start, statemod = statemod)

# compute initial likelihood to see if start is reasonable
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1, j = 1)
obj$get_par("sigma", k = 1, j = 1, s = 1:4)
obj$get_par("D")
obj$state()$delta() 

# predict membership for each individual 
pr <- obj$pr_state() 

# compare to truth
predstate <- matrix(0, nr = scrdat$n(), nc = 4)
for (i in 1:scrdat$n()) {
  predstate[i,] <- colSums(pr[[i]][,,10])
}
print(cbind(round(predstate, 2), cmix))




