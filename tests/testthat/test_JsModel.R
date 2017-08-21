context("Test JSModel object")

data("openpop_example_data")
time <- c(1, 2, 3, 5, 6, 9, 10, 11, 12, 13)
scr_dat <- ScrData$new(openpop_example_data$capture.history, 
                       openpop_example_data$mesh, time)

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            beta ~ 1, 
            phi ~ 1)

start <- list(lambda0 = 2.0, 
              sigma = 100, 
              beta = 0.2, 
              phi = 0.5, 
              D = 100)

mod <- JsModel$new(par, scr_dat, start, num_cores = 3)

test_that("data is loaded in model is equal to ScrObject", 
          expect_identical(mod$data(), scr_dat))

test_that("arrival probabilities sum to one", 
          expect_equal(sum(mod$get_par("beta")), 1))

test_that("pdet is computed correctly", 
          expect_less_than(mod$calc_pdet() -  0.739126, 1e-6))

test_that("initial distribution is correct", 
          expect_identical(mod$calc_initial_distribution() * scr_dat$n_meshpts(), 
                           matrix(c(1 - start$beta, start$beta, 0), nrow = scr_dat$n_meshpts(), ncol = 3, byrow = TRUE)))

b <- c(start$beta, (1 - start$beta) * (diff(time) / sum(diff(time))))
denom <- 1 - start$beta
pr_e <- rep(0, 9)
for (i in 1:9) {
  pr_e[i] <- b[i + 1] / denom 
  denom <- denom * (1 - pr_e[i])
}

test_that("pr_entry is correct", 
          expect_equal(mod$calc_pr_entry(), pr_e))

