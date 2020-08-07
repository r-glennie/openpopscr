#### Basic SCR example with pairwise non-euclidean distance
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 1)

# simulate data -----------------------------------------------------------
# set truth 
true_par <- list(D = 1000, lambda0 = 0.2, sigma = 20, noneuc = 0.5)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set non-euclidean pairwise distance effect
euc <- matrix(0, nrow = nrow(detectors), nc = nrow(mesh))
euc[1:10, 1:1500] <- 1 
sigma_mesh <- exp(log(true_par$sigma) + true_par$noneuc * euc)
plot(mesh[,1], mesh[,2], col = sigma_mesh[1,], pch = 19)
plot(mesh[,1], mesh[,2], col = sigma_mesh[20,], pch = 19)

# set number of occasions to simulate
n_occasions <- 5 

# simulator ---------------------------------------------------------------

# edit simulator to include/exclude covariates
simulate_survey <- function() {
  rownames(detectors) <- 1:nrow(detectors)
  D <- true_par$D
  lambda0 <- true_par$lambda0
  sigma <- true_par$sigma
  K <- n_occasions
  A <- nrow(mesh) * attr(mesh, "area") / 100
  N <- rpois(1, D*A)
  pt <- sample(1:nrow(mesh), size = N, replace = TRUE)
  x <- mesh[pt, 1]
  y <- mesh[pt, 2]
  
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
      er <- lambda0 * exp(-d2 / (2 * sigma_mesh[,pt[i]]^2))
      c <- rpois(length(er), er)
      if (any(c > 0)) {
        if (!seen[i]) {
          id[i] <- max(id) + 1
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
  
  scrdat <- ScrData$new(ch, mesh = mesh)
  return(scrdat)
}


# simulate data -----------------------------------------------------------

scrdat <- simulate_survey()

# openpopscr fit ----------------------------------------------------------

# create custom hazard function
hfn <- function(x, par) {
  sigma <- exp(log(par[[2]]) + euc * par[[3]])
  par[[1]]*exp(-(x^2)/(2*sigma^2))
}

# create detection function object 
noneuc_detfn <- DetFn$new(c("lambda0", "sigma", "noneuc"), 
                    fn = hfn, 
                    prob = FALSE, 
                    link2response = list("exp", "exp", "identity"), 
                    response2link = list("log", "log", "identity"))

# create formulae 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             noneuc ~ 1, 
             D ~ 1)

# set starting values  
start <- list(lambda0 = 0.2, sigma = 20, noneuc = 0, D = 1000)

# create the model object 
obj <- ScrModel$new(form, scrdat, start, detectfn = noneuc_detfn)

# compute initial likelihood 
obj$calc_llk()

# fit model
obj$fit()

# see model results 
obj

# get parameters on natural scale 
obj$get_par("lambda0", k = 1, j = 1)
obj$get_par("sigma", k = 1, j = 1)
obj$get_par("noneuc", k = 1, j = 1)
obj$get_par("D")
