#### Simple script to simulate SCR data 

## Set base parameters 
D <- 1000
lambda0 <- 1.0
sigma <-  20

## Survey setup 
# number of occasions
K <- 5
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

## Create covariates
detcov <- rnorm(nrow(detectors)) 
detcov <- detcov - mean(detcov)
lambda0 <- exp(log(lambda0) + 0.3 * detcov)

## Simulate activity centres
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
    er <- lambda0 * exp(-d2 / (2 * sigma^2))
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
scrdat$add_covariate("detcov", detcov, "j")

library(secr)
covariates(detectors) <- data.frame(detcov = detcov)
ch2 <- make.capthist(cap, detectors)
secrfit <- secr.fit(ch2, model = list(D ~ 1, lambda0 ~ detcov, sigma ~ 1), mask = mesh, detectfn="HHN")

#secrfit0 <- secr.fit(ch, model = list(D ~ 1, lambda0 ~ 1, sigma ~ 1), mask = mesh, detectfn="HHN")

