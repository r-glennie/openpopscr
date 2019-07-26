#### SCR with covariates on detection simulation 
library(openpopscr)
# set number of threads for parallel processing 
RcppParallel::setThreadOptions(numThreads = 30)

# setup simulations -------------------------------------------------------
set.seed(5392)
nsims <- 100 
ests <- vector(mode = "list", length = nsims)

# set truth 
true_par <- list(D = 1000, lambda0 = 0.2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# Create covariates
detcov <- rnorm(nrow(detectors)) 
detcov <- detcov - mean(detcov)
detcovf <- factor(sample(1:2, size = nrow(detectors), replace = TRUE))
timecov <- rnorm(n_occasions)
timecov <- timecov - mean(timecov)
timecovf <- factor(sample(1:3, size = n_occasions, replace = TRUE))
dettime <- matrix(rnorm(n_occasions*nrow(detectors)), nr = n_occasions, nc = nrow(detectors))

# create formulae for model to test  
form <- list(lambda0 ~ detcov, 
             sigma  ~ 1, 
             D ~ 1)


# simulator ---------------------------------------------------------------

# edit simulator to include/exclude covariates
simulate_survey <- function() {
  # set formulae here 
  lambda0 <- exp(log(true_par$lambda0) + 0.2 * detcov)
  sigma <- exp(log(true_par$sigma))
  
  # set ER function here 
  compute_er <- function(d2) {lambda0 * exp(-d2/(2*sigma^2))}
  
  # Simulate activity centres
  A <- nrow(mesh) * attr(mesh, "area") / 100
  N <- rpois(1, true_par$D*A)
  pt <- sample(1:nrow(mesh), size = N, replace = TRUE)
  x <- mesh[pt, 1]
  y <- mesh[pt, 2]
  
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
      er <- compute_er(d2)
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
  scrdat$add_covariate("detcov", detcov, "j")
  scrdat$add_covariate("detcovf", detcovf, "j")
  scrdat$add_covariate("timecov" , timecov, "k")
  scrdat$add_covariate("timecovf" , timecovf, "k")
  scrdat$add_covariate("dettime", dettime, "kj")
  return(scrdat)
}

# simulation --------------------------------------------------------------
progbar <- utils::txtProgressBar(min = 0, max = nsims, style = 3)  

for (sim in 1:nsims) {
  
  ## progress bar
  Sys.sleep(0.1) 
  utils::setTxtProgressBar(progbar, sim)
  
  # simulate ScrData 
  scrdat <- simulate_survey()
  
  # get starting values for numerical optimiser  
  start <- get_start_values(scrdat)
  
  # create the model object 
  obj <- ScrModel$new(form, scrdat, start, print = FALSE)
  
  # fit model
  obj$fit()
  
  # store results
  ests[[sim]] <- obj$estimates()$par
  
}

mu <- sapply(ests, FUN = function(x){x[,1]})
lcl <- sapply(ests, FUN = function(x){x[,3]})
ucl <- sapply(ests, FUN = function(x){x[,4]})
# distributions 
summary(t(mu))
# confidence interval coverage 
sum(lcl[1,] < log(true_par$lambda0) & ucl[1,] > log(true_par$lambda0))
sum(lcl[2,] < 0.2 & ucl[2,] > 0.2)
sum(lcl[3,] < log(true_par$sigma) & ucl[3,] > log(true_par$sigma), na.rm = TRUE)
sum(lcl[4,] < log(true_par$D) & ucl[4,] > log(true_par$D))
