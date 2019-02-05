# SCR example 
library(openpopscr)
library(secr)
RcppParallel::setThreadOptions(numThreads = 4)

# set seed 
set.seed(58229)

# set truth ---------------------------------------------------------------

# parameters
true_par <- list(D = 1000, lambda0 = 0.5, sigma = 20)

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# set number of occasions to simulate
n_occasions <- 5 


# simulation setup --------------------------------------------------------

# number of simulations 
nsims <- 10 

# save openpopscr models 
obj_list <- vector(mode = "list", length = nsims)

form <- list(lambda0 ~ 1, 
             sigma  ~ 1)

# run simulation ----------------------------------------------------------

for (sim in 1:nsims) {

  # simulate ScrData 
  scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh)

  #  get starting values 
  start <- get_start_values(scrdat)

  # create model object  
  obj <- ScrModel$new(form, scrdat, start, print = FALSE)

  # fit model
  obj$fit()

  # save model 
  obj_list[[sim]] <- obj 
  
  # save results to file 
  saveRDS(obj_list, "scr_basic_results.Rds")

}
