context("JsModel test")

# set truth 
true_par <- list(D = 500, lambda0 = 0.2, sigma = 30, phi = 0.5, beta = 0.3)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5 
# set N 
N <- 100
# simulate ScrData 
scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, print = FALSE, seed = 18953)
# formula 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             beta ~ 1, 
             phi ~ 1, 
             D ~ 1)
# get automatic start values 
start <- get_start_values(scrdat, model = "JsModel")
# create object
obj <- JsModel$new(form, scrdat, start, print = FALSE)

test_that("Computes correct pdet", {
  expect_equal(signif(obj$calc_pdet(), 4), 0.4633)
})  

test_that("Initial distribution valid", {
  ini <- obj$calc_initial_distribution()
  expect_equal(ncol(ini), 3)
  expect_equal(sum(ini), start$D*scrdat$cell_area())
  expect_equal(sum(ini[,2]), start$beta*start$D*scrdat$cell_area())
  expect_equal(sum(ini[,3]), 0)
})

test_that("Tpms are valid", {
  tpms <- obj$calc_tpms()[-5]
  sums <- sapply(tpms, FUN = function(x){rowSums(x)})
  expect_true(all(sums == 1))
  neg <- sapply(tpms, FUN = function(x){any(x < -1e-10)})
  expect_true(all(neg == FALSE))
  pos <- sapply(tpms, FUN = function(x){any(x > 1 + 1e-10)})
  expect_true(all(pos == FALSE))
})

test_that("Llk has correct valid", {
  expect_equal(signif(obj$calc_llk(), 4), -400.9)
})

test_that("Models fitting works", {
  suppressWarnings(obj$fit())
  expect_equal(signif(obj$estimates()$par, 4), c(-1.867,3.385,0.371,-0.6073,6.352,0.2782,0.131,0.5229,0.5659,0.2755,-2.412,3.128,-0.6539,-1.716,5.812,-1.322,3.642,1.396,0.5019,6.892), check.attributes = FALSE)
})

n_occasions2 <- 10

# set primary periods 
primary <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2))

# simulate ScrData 
scrdat2 <- simulate_js_openscr(true_par, 
                               n_occasions2, 
                               detectors, 
                               mesh, 
                               primary = primary, 
                               print = FALSE,
                               seed = 19482)

# get starting values 
start2 <- get_start_values(scrdat, model = "JsModel")

# create model object 
oo <- JsModel$new(form, scrdat2, start2, print = FALSE)
suppressWarnings(oo$fit())

test_that("Robust design model works", {
  expect_equal(signif(oo$estimates()$par, 4), c(-1.483,3.225,0.393,-1.051,6.862,0.1558,0.08106,0.4481,0.4968,0.2817,-1.789,3.066,-0.4852,-2.025,6.31,-1.178,3.384,1.271,-0.07764,7.414), check.attributes = FALSE)
})
