#### SCR with mixture on detection simulation 
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 30)

# setup simulations -------------------------------------------------------
set.seed(15583)
nsims <- 100
ests <- vector(mode = "list", length = nsims)

# set truth
D <- 1000
lambda0 <- c(0.5, 0.5)
sigma <- c(20, 40)
delta <- c(0.3, 0.7)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# create formulae for model to test  
form <- list(lambda0 ~ 1, 
             sigma  ~ state, 
             D ~ 1)


# simulator ---------------------------------------------------------------

# edit simulator to include/exclude covariates
simulate_survey <- function() {
  # set ER function here 
  compute_er <- function(d2) {lambda0 * exp(-d2/(2*sigma^2))}
  
  # Simulate activity centres
  A <- nrow(mesh) * attr(mesh, "area") / 100
  N <- rpois(1, D*A)
  pt <- sample(1:nrow(mesh), size = N, replace = TRUE)
  x <- mesh[pt, 1]
  y <- mesh[pt, 2]
  
  ## Simulate mixture 
  mix <- sample(1:2, size = N, replace = TRUE, prob = delta)
  
  # Simulate survey 
  cap <- data.frame(session = numeric(), 
                    ID = numeric(), 
                    occasion = numeric(), 
                    trap = numeric())
  
  seen <- rep(FALSE, N)
  id <- rep(0, N)
  for (k in 1:n_occasions) {
    for (i in 1:N) {
      d2 <- (x[i] - detectors[,1])^2 + (y[i] - detectors[,2])^2
      er <- lambda0[mix[i]] * exp(-d2 / (2 * sigma[mix[i]]^2))
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
  if (max(cap$occasion) != n_occasions) cap <- rbind(cap, data.frame(session = 1, ID = "NONE", occasion = n_occasions, trap = 1))
  ch <- make.capthist(cap, detectors)
  scrdat <- ScrData$new(ch, mesh = mesh)
  return(scrdat)
}

# simulation --------------------------------------------------------------
progbar <- utils::txtProgressBar(min = 0, max = nsims, style = 3)  

for (sim in 1:nsims) {
  
  ## progress bar
  Sys.sleep(0.1) 
  utils::setTxtProgressBar(progbar, sim)
  
  test <-  try({
  # simulate ScrData 
  scrdat <- simulate_survey()
  
  # get starting values for numerical optimiser  
  start <- get_start_values(scrdat)
  
  # create state model
  statemod <- StateModel$new(data = scrdat, 
                             names = c("M", "F"), 
                             structure = matrix(c(".", "0", 
                                                  "0", "."), nr = 2, nc = 2, byrow = T), 
                             start = list(delta = c(0.3, 0.7), tpm = diag(2)))
  
  # create the model object 
  obj <- ScrModel$new(form, scrdat, start, statemod = statemod, print = FALSE)
  
  # fit model
  obj$fit()
  
  est <- list(est = obj$estimates()$par, delta = obj$state()$estimates())
  
  })
  
  if ("try-error" %in% class(test)) {
    ests[[sim]] <- NA
  } else {
    ests[[sim]] <- est 
  }
}

## Estimates 
mu <- sapply(ests, FUN = function(x){x$est[,1]})
lcl <- sapply(ests, FUN = function(x){x$est[,3]})
ucl <- sapply(ests, FUN = function(x){x$est[,4]})
# distributions 
summary(t(mu))
# confidence interval coverage 
sum(lcl[1,] < log(lambda0[1]) & ucl[1,] > log(lambda0[1]))
sum(lcl[2,] < log(sigma[1]) & ucl[2,] > log(sigma[1]))
sum(lcl[3,] < log(sigma[2]) - log(sigma[1]) & ucl[3,] > log(sigma[2]) - log(sigma[1]))
sum(lcl[4,] < log(D) & ucl[4,] > log(D))

## Mixture 
pmix <- sapply(ests, FUN = function(x){x$delta[1]})
plcl <- sapply(ests, FUN = function(x){x$delta[3]})
pucl <- sapply(ests, FUN = function(x){x$delta[3]})
summary(pmix)
sum(plcl < qlogis(delta[1]) & qlogis(delta[1]) < pucl)




