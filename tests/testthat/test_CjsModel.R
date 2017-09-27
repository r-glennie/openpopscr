context("Test CjsModel object")

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

mod <- CjsModel$new(par, scr_dat, start, num_cores = 3)

test_that("data is loaded in model is equal to ScrObject", 
          expect_identical(mod$data(), scr_dat))

test_that("initial distribution is correct", 
          expect_identical(mod$calc_initial_distribution() * scr_dat$n_meshpts(), 
                           matrix(c(1, 0), nrow = scr_dat$n_meshpts(), ncol = 2, byrow = TRUE)))



