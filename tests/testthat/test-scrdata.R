context("ScrData test")
library(secr)
set.seed(53838)
detectors <- make.grid(nx = 3, ny = 3, spacing = 20, detector = "multi")
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
popn <- sim.popn(D = 100, core = mesh, Ndist = "poisson", buffertype = "rect")
ch <- sim.capthist(detectors,  
             popn = popn, 
             detectfn = "HHN", 
             detectpar = list(lambda0 = 0.5, 
                              sigma = 20), 
             noccasions = 5, 
             renumber = FALSE)

test_that("ScrData checks input", {
  expect_error(ScrData$new(as.vector(ch), mesh), "Invalid capture history object.")
  expect_error(ScrData$new(ch, as.matrix(mesh)), "Invalid mesh object.")
  expect_error(ScrData$new(ch, mesh, time = rep("a", 5)), "Time is not a numeric vector.")
  expect_error(ScrData$new(ch, mesh, time = rep(1, 4)), "Length of time vector not equal to number of occasions.")
  expect_error(ScrData$new(ch, mesh, time = rep(1, 10)), "Length of time vector not equal to number of occasions.")      
  expect_error(ScrData$new(ch, mesh, time = c(1,2,3,5,4)), "Time must be an increasing vector of numbers.") 
  expect_error(ScrData$new(ch, mesh, primary = rep("a", 5)), "Primary is not a numeric vector.")
  expect_error(ScrData$new(ch, mesh, primary = rep(1, 4)), "Length of primary vector not equal to number of occasions.")
  expect_error(ScrData$new(ch, mesh, primary = rep(1, 10)), "Length of primary vector not equal to number of occasions.")      
  expect_error(ScrData$new(ch, mesh, primary = c(1,2,3,5,4)), "Primary must be an increasing vector of numbers.") 
  expect_error(ScrData$new(ch, mesh, primary = c(1,2,4,5,5)), "Primary must contain integers 1:maximum primary.") 
  expect_error(ScrData$new(ch, mesh, primary = c(1,2,4,5,6)), "Primary must contain integers 1:maximum primary.")
  expect_error(ScrData$new(ch, mesh, primary = c(1,2,4.2,5,6)), "Primary must be integer labels.")
  }
) 

test_that("ScrData recognises detector types", {
  attributes(ch)$traps <- make.grid(nx = 3, ny = 3, spacing = 20, detector = "count")
  expect_equal(ScrData$new(ch, mesh)$detector_type(), 1)
  attributes(ch)$traps <- make.grid(nx = 3, ny = 3, spacing = 20, detector = "proximity")
  expect_equal(ScrData$new(ch, mesh)$detector_type(), 2)
  attributes(ch)$traps <- make.grid(nx = 3, ny = 3, spacing = 20, detector = "multi")
  expect_equal(ScrData$new(ch, mesh)$detector_type(), 3)
  attributes(ch)$traps <- make.grid(nx = 3, ny = 3, spacing = 20, detector = "single")
  suppressWarnings(expect_equal(ScrData$new(ch, mesh)$detector_type(), 3))
  expect_warning(ScrData$new(ch, mesh)$detector_type(), "Single-catch detectors treated as multi-catch detectors.")
}
)
          
test_that("ScrData retains trap usage", {
  attributes(detectors)$usage <- matrix(1:45, nr = 9, nc = 5)
  attributes(ch)$traps <- detectors
  expect_equal(usage(ScrData$new(ch, mesh)$traps()), usage(detectors))
}
)         
 
test_that("Scrdata returns correct times", {
  expect_equal(ScrData$new(ch, mesh)$time(), 1:5)
  expect_equal(ScrData$new(ch, mesh, time = 1:5)$time(), 1:5)
  expect_equal(ScrData$new(ch, mesh, time = c(1, 1.5, 2.4, 3.55, 6))$time(), c(1, 1.5, 2.4, 3.55, 6))
  expect_equal(ScrData$new(ch, mesh, primary = 1:5)$primary(), 1:5)
  expect_equal(ScrData$new(ch, mesh, primary = c(1,1,2,2,3))$primary(), c(1,1,2,2,3))
  expect_equal(ScrData$new(ch, mesh)$primary(), seq(1,5))
}         
)

test_that("ScrData computes correct statistics", {
  dat <- ScrData$new(ch, mesh)
  expect_equal(dat$n(), 138)
  expect_equal(dat$n_occasions(), 5)
  expect_equal(dat$n_traps(), 9)
  expect_equal(dat$n_meshpts(), 3496)
  expect_equal(dat$distances()[1,5], sqrt(0.625^2+98.125^2))
  expect_equal(dat$distances()[5,5], sqrt((20+0.625)^2+(20+98.125)^2))
  expect_equal(dat$area(), nrow(mesh) * attr(mesh, "area") / 100)
  expect_equal(dat$area(), nrow(mesh) * dat$cell_area())
  expect_equal(dat$encrange(), RPSV(ch))
  dat <- ScrData$new(ch, mesh, primary = c(1,1,1,2,3))
  expect_equal(dat$n_occasions(), 3)
  expect_equal(dat$n_primary(), 3)
  expect_equal(dat$n_occasions("all"), 5)
  expect_equal(dat$n_secondary(), c(3, 1, 1))
}
)

test_that("ScrData checks covariates", {
  dat <- ScrData$new(ch, mesh)
  expect_error(dat$add_covariate(1, rep(1:9), "j"), "Covariate name must be a character string.")
  expect_error(dat$add_covariate("t", rep(1:9), "j"), "Covariate with that name already exists.")
  expect_error(dat$add_covariate("cov1", rep(1:9), "j2"), "Invalid covariate type.")
  expect_error(dat$add_covariate("cov2", rep("a", 5), "j"), "Invalid covariate, must be factor or numeric.")
  expect_error(dat$add_covariate("cov3", 1:4, "j"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov4", matrix(0, nr = 9, nc = 9), "j"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov6", 1:6, "k"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov7", matrix(0, nr = 3, nc = 3), "k"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov8", 1:45, "kj"))
  expect_error(dat$add_covariate("cov9", matrix(0, nr = 9, nc = 5), "kj"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov10", 1:45, "km"))
  expect_error(dat$add_covariate("cov11", matrix(0, nr = 3496, nc = 5), "km"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov12", matrix(0, nr = 5 , nc = 3495), "km"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov13", 1:3494, "m"), "Invalid covariate.")
  expect_error(dat$add_covariate("cov14", matrix(0, nr = 59, nc = 59), "j"), "Invalid covariate.")
})

test_that("ScrData loads covariates", {
  dat <- ScrData$new(ch, mesh)
  dat$add_covariate("detcov", 1:9, "j")
  dat$add_covariate("timecov", 1:5, "k")
  dat$add_covariate("dettime", matrix(1:45, nr = 5, nc = 9), "kj")
  dat$add_covariate("timemesh", matrix(1:(5*3496), nr = 5, nc = 3496), "km")
  dat$add_covariate("mesh", 1:3496, "m")
  expect_equal(names(dat$covs()), c("t", "T", "primary", "Primary", "x", "y", "detcov", "timecov", "dettime", "timemesh", "mesh"))
  expect_equal(dat$covs(j = 1)$detcov, 1)
  expect_equal(dat$covs(j = 4)$detcov, 4)
  expect_equal(dat$covs(k = 1)$timecov, 1)
  expect_equal(dat$covs(k = 5)$timecov, 5)
  expect_equal(dat$covs(k = 1, j = 2)$dettime, 6)
  expect_equal(dat$covs(k = 3, j = 4)$dettime, 18)
  expect_equal(dat$covs(k = 1, m = 2)$timemesh, 6)
  expect_equal(dat$covs(k = 3, m = 4)$timemesh, 18)
  expect_equal(dat$covs(m = 5)$mesh, 5)
  expect_equal(dat$covs(k = c(1, 5))$timecov, c(1, 5))
  expect_equal(dat$covs(j = c(1, 9))$detcov, c(1, 9))
  expect_equal(dat$covs(j = c(3,4,6), k = c(1, 2, 3))$dettime, matrix(c(11, 12, 13, 16, 17, 18, 26, 27, 28), nr = 3))
}
)

