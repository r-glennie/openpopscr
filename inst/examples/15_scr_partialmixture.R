#### Basic partially observed mixture SCR example
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------
set.seed(53919)

## Set parameters 
D <- 1000
lambda0 <- c(0.5, 0.5)
sigma <-  c(20, 40)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
K <- 5

## Simulate activity centres
A <- nrow(mesh) * attr(mesh, "area") / 100
N <- rpois(1, D*A)
pt <- sample(1:nrow(mesh), size = N, replace = TRUE)
x <- mesh[pt, 1]
y <- mesh[pt, 2]

## Simulate mixture 
nstates <- 2 
delta <- c(0.7, 0.3)
mix <- sample(1:2, size = N, replace = TRUE, prob = delta)
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
si <- rep(NA, scrdat$n())
si[1:20] <- cmix[1:20]
# add to data 
scrdat$add_covariate("sex", si, "i")

# openpopscr fit ----------------------------------------------------------

## create state model 
statemod <- StateModel$new(data = scrdat, 
                           names = c("female", "male"), 
                           structure = matrix(c(".", "0", 
                                                "0", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), tpm = diag(2)), 
                           cov = data.frame(sex = factor(c(1, 2))))

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ sex,
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
obj$get_par("sigma", k = 1, j = 1, s = 1)
obj$get_par("sigma", k = 1, j = 1, s = 2)
obj$get_par("D")
statemod$delta()
