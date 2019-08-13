context("CjsModel test")

# set truth 
true_par <- list(lambda0 = 0.2, sigma = 30, phi = 0.7)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5 
# set N 
N <- 100
# simulate ScrData 
scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, print = FALSE, seed = 19295)
# formula 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             phi ~ 1)
# get automatic start values 
start <- get_start_values(scrdat, model = "CjsModel")
# create object
obj <- CjsModel$new(form, scrdat, start, print = FALSE)

test_that("Computes correct inipdet", {
  prcap <- obj$calc_pr_capture()
  inipdet <- obj$calc_initial_pdet(prcap)
  expect_equal(signif(inipdet, 4), c(7.978e-06,0.6002,0.0002516,0.6002,0.6002,0.0007688,0.001527,0.6002,0.6002,0.6002,0.6002,0.0002516,0.6002,9.448e-05))
})  

test_that("Initial distribution valid", {
  ini <- obj$calc_initial_distribution()
  expect_equal(ncol(ini), 2)
  expect_equal(sum(ini), 1)
  expect_equal(sum(ini[,2]), 0)
})

test_that("Tpms are valid", {
  tpms <- obj$calc_tpms()[-5]
  sums <- sapply(tpms, FUN = function(x){rowSums(x)})
  expect_true(all(sums == 1))
  neg <- sapply(tpms, FUN = function(x){any(x < 0)})
  expect_true(all(neg == FALSE))
  pos <- sapply(tpms, FUN = function(x){any(x > 1)})
  expect_true(all(pos == FALSE))
})

test_that("Llk has correct valid", {
  expect_equal(signif(obj$calc_llk(), 4), -124.5)
})

test_that("Entry occasions are correct", {
  ex <- obj$entry()
  names(ex) <- NULL
  expect_equal(ex, c(0,1,0,1,2,0,0,3,1,3,1,0,4,0))
})

test_that("Models fitting works", {
  obj$fit()
  expect_equal(signif(obj$estimates()$par, 4), c(-0.4232,3.838,1.129,0.6615,0.1132,0.5485,-1.72,3.616,0.05399,0.8733,4.059,2.204), check.attributes = FALSE)
})

n_occasions2 <- 10

# set primary periods 
primary <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2))

# set N 
N2 <- 200

# simulate ScrData 
scrdat2 <- simulate_cjs_openscr(true_par, 
                               N2, 
                               n_occasions2, 
                               detectors, 
                               mesh, 
                               primary = primary, 
                               print = FALSE,
                               seed = 19482)

# get starting values 
start2 <- get_start_values(scrdat, model = "CjsModel")

# create model object 
oo <- CjsModel$new(form, scrdat2, start2, print = FALSE)
oo$fit()

test_that("Robust design model works", {
  expect_equal(signif(oo$estimates()$par, 4), c(-1.624,3.338,0.7767,0.09872,0.05648,0.2278,-1.818,3.228,0.3303,-1.431,3.449,1.223), check.attributes = FALSE)
})
