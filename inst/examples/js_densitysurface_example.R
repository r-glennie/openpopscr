# Jolly-Seber SCR density surface example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 6)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20, phi = 0.8, beta = 0.3, sd = 10)

# make detectors array 
detectors <- secr::make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- secr::make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# spatial density effect
x.std <- scale(mesh[,1])
y.std <- scale(mesh[,2])
logD <- log(true_par$D) + 0.5 * x.std + 0.8 * y.std
ihp <- exp(logD) / true_par$D

# set number of occasions to simulate
n_occasions <- 10

# simulate ScrData 
scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, ihp = ihp, move = TRUE)

# openpopscr fit ----------------------------------------------------------

form <- list(lambda0 ~ 1, 
            sigma ~ 1,
            beta ~ 1, 
            phi ~ 1,
            D ~ 1)

start <- get_start_values(scrdat, model = "JsModel")

obj <- JsModel$new(form, scrdat, start)

# fit model
obj$fit()

# plot density surface
D <- obj$get_par("D", m = 1:scrdat$n_meshpts()) * scrdat$cell_area()
scrdat$plot_mesh(D)

form2 <- list(lambda0 ~ 1, 
             sigma ~ 1,
             beta ~ 1, 
             phi ~ 1,
             D ~ x)

obj2 <- JsModel$new(form2, scrdat, start)

# fit model
obj2$fit()

# plot density surface
D <- obj2$get_par("D", m = 1:scrdat$n_meshpts()) * scrdat$cell_area()
scrdat$plot_mesh(D)


fit <- openCR.fit(scrdat$capthist(), mask = scrdat$mesh(), type = "secrD", trace = TRUE)

# see model results 
obj

# get parameters on natural scale 
obj2$get_par("lambda0", k = 1)
obj2$get_par("sigma", k = 1)
obj2$get_par("phi", k = 1)
obj2$get_par("beta", k = 1)
obj2$get_par("D")


form3 <- list(lambda0 ~ 1, 
              sigma ~ 1,
              beta ~ 1, 
              phi ~ 1,
              sd ~ 1, 
              D ~ x + y)
start3 <- get_start_values(scrdat, model = "JsTransientModel")
obj3 <- JsTransientModel$new(form3, scrdat, start3)
obj3$calc_llk()
obj3$fit()


