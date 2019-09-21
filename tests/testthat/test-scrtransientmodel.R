context("ScrTransientModel test")

# set truth 
true_par <- list(D = 1000, lambda0 = 0.5, sigma = 20, sd = 10)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5 
# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, move = TRUE, print = FALSE, seed = 15483)
# formula 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             sd ~ 1, 
             D ~ 1)
# get automatic start values 
start <- get_start_values(scrdat, model = "ScrTransientModel")

test_that("Basic model computes correct llk", {
  obj <- ScrTransientModel$new(form, scrdat, start, print = FALSE)
  expect_equal(signif(obj$calc_llk(), 6), -1292.12)
}  
)

