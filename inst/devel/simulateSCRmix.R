#### Simple script to simulate SCR data 
## Set parameters 
D <- 500
lambda0 <- c(0.5, 0.5)
sigma <-  c(50, 20)

## Survey setup 
# number of occasions
K <- 5
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
rownames(detectors) <- 1:nrow(detectors)
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

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
#tpm <- matrix(c(0.8, 0.2, 
#                0.3, 0.7), nr = nstates, nc = nstates, byrow = TRUE)
tpm <- diag(2)

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
    mix[i] <- sample(1:nstates, size = 1, prob = tpm[mix[i],])
  }
}
if (max(cap$occasion) != K) cap <- rbind(cap, data.frame(session = 1, ID = "NONE", occasion = K, trap = 1))
ch <- make.capthist(cap, detectors)

scrdat <- ScrData$new(ch, mesh = mesh)



## create state model 
statemod <- StateModel$new(data = scrdat, 
                           names = c("male", "female"), 
                           structure = matrix(c(".", "0", 
                                                "0", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), tpm = diag(2)))#tpm = matrix(c(0.8, 0.2,
                                                                            #0.2, 0.8), nr = 2, nc = 2, byrow = T)))
                           #delta_fixed = c(TRUE, TRUE))

form <- list(lambda0 ~ 1, 
             sigma ~ state, 
             D ~ 1)

start <- get_start_values(scrdat)

#mod <- ScrModel$new(form, scrdat, start, statemod = statemod)

#mod$calc_llk()

#mod$fit()

## some observed sexes 
s <- ifelse(cmix == 1, "male", "female")
si <- vector(mode = "list", length = scrdat$n())
for (i in 1:20) {
  si[[i]] <- s[i]
}

scrdat$add_covariate("state", si, "i")

mod <- ScrModel$new(form, scrdat, start, statemod = statemod)

mod$calc_llk()

mod$fit()


library(secr)
ch <- scrdat$capthist()
covariates(ch) <- data.frame(sex = factor(c(s[1:20], rep(NA, length(s) - 20))))
secrfit <- secr.fit(ch, 
                    mask = scrdat$mesh(), 
                    detectfn = "HHN", 
                    model = list(sigma ~ h2), 
                    hcov = "sex")



