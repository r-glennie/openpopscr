#### Basic mixture SCR example
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------
set.seed(53919)

## Set parameters 
D <- 1000
lambda0 <- c(1.0, 1.0)
sigma <-  c(20, 40)

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
nstates <- 2 
delta <- c(0.7, 0.3)
mix <- sample(1:2, size = N, replace = TRUE, prob = delta)

## Simulate survey 
cap <- data.frame(session = numeric(), 
                  ID = numeric(), 
                  occasion = numeric(), 
                  trap = numeric())

seen <- rep(FALSE, N)
id <- rep(0, N)
obsmix <- NULL
for (k in 1:K) {
  for (i in 1:N) {
    d2 <- (x[i] - detectors[,1])^2 + (y[i] - detectors[,2])^2
    er <- lambda0[mix[i]] * exp(-d2 / (2 * sigma[mix[i]]^2))
    c <- rpois(length(er), er)
    if (any(c > 0)) {
      if (!seen[i]) {
        id[i] <- max(id) + 1
        seen[i] <- TRUE
        obsmix <- c(obsmix, mix[i])
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

# openpopscr fit ----------------------------------------------------------

## create state model 
statemod <- StateModel$new(data = scrdat, 
                           names = c("A", "B"), 
                           structure = matrix(c(".", "0", 
                                                "0", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), tpm = diag(2)))

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ state,
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
obj$state()$delta()

# predict state 
pr <- obj$pr_state()
predstate <- matrix(0, nr = scrdat$n(), nc = 2)
for (i in 1:scrdat$n()) {
  predstate[i,] <- colSums(pr[[i]][,,10])
}
print(cbind(round(predstate, 2), obsmix))

predmix <- rep(0, scrdat$n())
for (i in 1:scrdat$n()) predmix[i] <- which.max(predstate[i,])
sum(predmix == obsmix) / scrdat$n()
