context("CjsTransientModel test")

# set truth 
true_par <- list(lambda0 = 2, sigma = 20, phi = 0.8, sd = 10)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5 
# number of individuals
N <- 100
# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, move = TRUE, print = FALSE, seed = 58348)
# formula 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             phi ~ 1, 
             sd ~ 1)
# get automatic start values 
start <- get_start_values(scrdat, model = "CjsTransientModel")

test_that("Basic model computes correct llk", {
  obj <- CjsTransientModel$new(form, scrdat, start, print = FALSE)
  expect_equal(signif(obj$calc_llk(), 6), -265.489)
}  
)

