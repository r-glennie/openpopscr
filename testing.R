# testing

library(openpopscr)
data("openpop_example_data")
impact <- factor(c(rep(0, 5), rep(1, 5)))

scr_dat <- ScrData$new(openpop_example_data$capture.history, 
                       openpop_example_data$mesh,
                       cov = list(impact = impact), 
                       cov_type = c("k"))


par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            beta ~ 1, 
            phi ~ impact)

start <- list(lambda0 = 2.0, 
              sigma = 100, 
              beta = 0.2, 
              phi = 0.5, 
              D = 100)

mod <- JsModel$new(par, scr_dat, start, num_cores = 3)

mod$fit()

saveRDS(mod, "fitted_mod_impact.Rds")
