context("Test SCR Data object")

data("openpop_example_data")
time <- c(1, 2, 3, 5, 6, 9, 10, 11, 12, 13)
covs <- list(x = 2, 
             xj = runif(49),
             xk = runif(10), 
             xm = runif(3760), 
             xjk = matrix(runif(49 * 10), nr = 49, nc = 10), 
             xjm = matrix(runif(49 * 3760), nr = 49, nc = 3760), 
             xkm = matrix(runif(10 * 3760), nr = 10, nc = 3760), 
             xjkm = array(runif(10 * 40 * 3760), dim = c(49, 10, 3760)))
cov_type <- c("1","j", "k", "m", "jk", "jm", "km", "jkm")
scr_dat <- ScrData$new(openpop_example_data$capture.history, 
                       openpop_example_data$mesh, 
                       time,
                       covs,
                       cov_type)
c <- scr_dat$covs(j = 2, k = 7, m = 1001)
new_cap <- openpop_example_data$capture.history
secr::usage(secr::traps(new_cap)) <- matrix(1, nr = 49, nc = 10)
dist <- scr_dat$distances()

test_that("capture history is equal to raw data", 
          expect_identical(scr_dat$capthist(), new_cap)
)

test_that("mesh is equal to raw mesh", 
          expect_identical(scr_dat$mesh(), openpop_example_data$mesh)
)

test_that("time vector is equal to raw time vector", 
          expect_identical(scr_dat$time(), time))

test_that("covariates are loaded correctly", { 
          expect_equal(c$xj, covs$xj[2])
          expect_equal(c$xk, covs$xk[7])
          expect_equal(c$xm, covs$xm[1001])
          expect_equal(c$xjk, covs$xjk[2, 7])
          expect_equal(c$xjm, covs$xjm[2, 1001])
          expect_equal(c$xkm, covs$xkm[7, 1001])
          expect_equal(c$xjkm, covs$xjkm[2, 7, 1001]) })

test_that("n = # of individuals is correct", 
          expect_equal(scr_dat$n(), 84))

test_that("n_occasions = # of occasions", 
          expect_equal(scr_dat$n_occasions(), 10))

test_that("n_traps = # of traps", 
          expect_equal(scr_dat$n_traps(), 49))

test_that("n_meshpts = # of mesh points", 
          expect_equal(scr_dat$n_meshpts(), 3884))

test_that("area = area of mesh in sq km", 
          expect_lt(scr_dat$area() - 1.365469, 1e-6))

test_that("check a few distance calculations", {
          expect_equal(round(dist[1, 50]), 276)
          expect_equal(round(dist[10, 37]), 566)})
