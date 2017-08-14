context("Test SCR Data object")

data("openpopscr_example_data")
time <- c(1, 2, 3, 5, 6, 9, 10, 11, 12, 13)
scr_dat <- ScrData$new(openpopscr_example_data$capture.history, 
                       openpopscr_example_data$mesh, time)
new_cap <- openpopscr_example_data$capture.history
secr::usage(secr::traps(new_cap)) <- matrix(1, nr = 49, nc = 10)
dist <- scr_dat$distances()

test_that("capture history is equal to raw data", 
          expect_identical(scr_dat$capthist(), new_cap)
)

test_that("mesh is equal to raw mesh", 
          expect_identical(scr_dat$mesh(), openpopscr_example_data$mesh)
)

test_that("time vector is equal to raw time vector", 
          expect_identical(scr_dat$time(), time))

test_that("n = # of individuals is correct", 
          expect_equal(scr_dat$n(), 75))

test_that("n_occasions = # of occasions", 
          expect_equal(scr_dat$n_occasions(), 10))

test_that("n_traps = # of traps", 
          expect_equal(scr_dat$n_traps(), 49))

test_that("n_meshpts = # of mesh points", 
          expect_equal(scr_dat$n_meshpts(), 3760))

test_that("area = area of mesh in sq km", 
          expect_equal(scr_dat$area(), 0.094))

test_that("check a few distance calculations", {
          expect_equal(dist[1, 50], sqrt(47.5^2 + 92.5^2))
  expect_equal(dist[10, 37], sqrt((40 + 17.5)^2 + (20 + 92.5)^2))
}

)