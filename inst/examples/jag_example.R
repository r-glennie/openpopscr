## Jaguar example
library(openpopscr)


# load data ---------------------------------------------------------------

jag_data <- readRDS("inst/examples/jag_dat.Rds")


# stationary model --------------------------------------------------------

stat_par <- list(lambda0 ~ t, 
            sigma ~ t, 
            beta ~ 1, 
            phi ~ 1)

stat_start <- list(lambda0 = 0.05, 
              sigma = 5000, 
              beta = 0.2, 
              phi = 0.8, 
              D = 0.05)


stat_mod <- JsModel$new(stat_par, jag_data, stat_start, num_cores = 4)

stat_mod$fit()

stat_mod

stat_mod$get_par("lambda0")
stat_mod$get_par("sigma")
stat_mod$get_par("phi", k = 1)
stat_mod$get_par("beta", k = 1)
stat_mod$get_par("D")

# transient model ---------------------------------------------------------

trans_par <- list(lambda0 ~ t, 
            sigma ~ t, 
            beta ~ 1, 
            phi ~ 1, 
            sd ~ 1)

trans_start <- list(lambda0 = 0.05, 
              sigma = 3000, 
              beta = 0.2, 
              phi = 0.8, 
              sd = 3000, 
              D = 0.05)


trans_mod <- JsTransientModel$new(trans_par, jag_data, trans_start, num_cores = 4, stepmax = 1)

trans_mod$calc_llk()

trans_mod$fit()

trans_mod

trans_mod$get_par("lambda0", k = 1)
trans_mod$get_par("sigma", k = 1)
trans_mod$get_par("phi", k = 1)
trans_mod$get_par("beta", k = 1)
trans_mod$get_par("sd", k = 1)
trans_mod$get_par("D")
