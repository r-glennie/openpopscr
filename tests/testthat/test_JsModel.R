context("Test JSModel object")

data("openpopscr_example_data")
time <- c(1, 2, 3, 5, 6, 9, 10, 11, 12, 13)
scr_dat <- ScrData$new(openpopscr_example_data$capture.history, 
                       openpopscr_example_data$mesh, time)
par <- list(lambda0 = 1, sigma = 10, beta = 0.2, phi = 0.8, D = 100)
mod <- JsModel$new(scr_dat, par, num_cores = 3)

test_that("data is loaded in model is equal to ScrObject", 
          expect_identical(mod$data(), scr_dat))

test_that("parameters is equal to initial value", 
          expect_identical(mod$par(), par))

test_that("pdet is computed correctly", 
          expect_equal(mod$calc_pdet(), 0.2373508))

test_that("initial distribution is correct", 
          expect_identical(mod$calc_initial_distribution() * scr_dat$n_meshpts(), 
                           matrix(c(1 - par$beta, par$beta, 0), nrow = scr_dat$n_meshpts(), ncol = 3, byrow = TRUE)))

b <- c(par$beta, (1 - par$beta) * rep(diff(time) / sum(diff(time))))
denom <- 1 - par$beta
pr_e <- rep(0, 9)
for (i in 1:9) {
  pr_e[i] <- b[i + 1] / denom 
  denom <- denom * (1 - pr_e[i])
}

test_that("pr_entry is correct", 
          expect_equal(mod$calc_pr_entry(), pr_e))

