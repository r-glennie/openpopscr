context("ScrModel test")

# set truth 
true_par <- list(D = 300, lambda0 = 0.5, sigma = 20)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5 
# simulate ScrData 
scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, seed = 25939, print = FALSE)

# model setup
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ 1)
start <- list(lamdba0 = 0.5, sigma = 20, D = 300)

test_that("ScrModel checks inputs", {
  expect_error(ScrModel$new(form, scrdat$capthist(), start, print = FALSE), "Must supply a ScrData object for data.")  
  expect_error(ScrModel$new(form[[1]], scrdat, start, print = FALSE), "Invalid list of formulae.")
  expect_error(ScrModel$new(list("lambda0 ~ 1", "sigma ~ 1", "D ~ 1"), scrdat, print = FALSE), "Invalid list of formulae.")
  expect_error(ScrModel$new(form, scrdat, c(0.5, 20, 300), print = FALSE), "start must be a list of starting values for each parameter.")
  expect_error(ScrModel$new(form, scrdat, start, print = "TRUE"), "print must be TRUE or FALSE.")
})

test_that("Basic model works", {
  mod <- suppressWarnings(ScrModel$new(form, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.7404,2.737,5.654,0.2864,0.1151,0.3412,-1.302,2.511,4.985,-0.179,2.962,6.323),
                                                      nr = 3, nc = 4), check.attributes = FALSE)
})

test_that("Detection covariates model works", {
  ## detcov 
  scrdat$add_covariate("detcov", c(-0.3625,-0.4869,1.678,0.2306,-1.673,-1.118,0.4265,
                                   -2.006,0.5957,-0.4685,0.823,-1.23,-0.2205,1.673,
                                   1.341,-1.483,0.2875,-1.081,1.014,-0.3805,1.369,
                                   -1.407,-0.4858,-0.1298,2.413,2.015,1.099,
                                   -1.879,-0.1079,2.272,-1.377,0.5986,-0.6217,
                                   0.2381,-0.6856,-1.245,-0.09166,0.03824,-0.7979,
                                   -0.2653,0.6487,-0.09794,-0.3511,-0.8161,0.5126,
                                   0.7343,-0.5492,-1.074,-0.2868), "j")
  formdetcov <- list(lambda0 ~ detcov, sigma ~ 1, D ~ 1)
  mod <- suppressWarnings(ScrModel$new(formdetcov, scrdat, start, print = FALSE))
  suppressWarnings(suppressWarnings(mod$fit()))
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.8915,0.2838,2.762,5.668,0.2925,0.1596,0.1173,0.3421,-1.465,-0.02902,2.532,4.997,-0.3183,0.5966,2.992,6.338),
                                                      nr = 4, nc = 4), check.attributes = FALSE)
  ## timecov
  scrdat$add_covariate("timecov", c(-0.7987,-0.1283,-1.542,0.248,-0.07588), "k")
  formtimecov <- list(lambda0 ~ 1, sigma ~ timecov, D ~ 1)
  mod <- suppressWarnings(ScrModel$new(formtimecov, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.7107,2.635,-0.1294,5.67,0.2848,0.1344,0.09602,0.3409,-1.269,2.371,-0.3176,5.002,-0.1525,2.898,0.05876,6.338),
                                                      nr = 4, nc = 4), check.attributes = FALSE)
  
  ## mixed detcov and timecov in formula 
  formdettimecov <- list(lambda0 ~ timecov + detcov, sigma ~ 1, D ~ 1)
  mod <- suppressWarnings(ScrModel$new(formdettimecov, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.9784,-0.1751,0.2838,2.762,5.668,0.3236,0.2597,0.1596,0.1173,0.3421,-1.613,-0.6841,-0.02902,2.532,4.997,-0.3441,0.334,0.5966,2.992,6.338),
                                                      nr = 5, nc = 4), check.attributes = FALSE)
  
  ## det and time cov 
  scrdat$add_covariate("dettime", matrix(c(-1.125,-0.3552,-0.7712,0.06758,1.097,0.4569,-0.9512,0.1436,1.354,0.9711,-0.4308,-0.1576,-0.3884,-0.0782,-0.9711,0.6355,0.7569,1.606,0.0005856,0.5376,0.9206,-0.0548,0.9233,0.9337,-0.8067,0.05913,-0.6932,1.272,-0.4794,0.0427,0.6546,0.02216,-1.147,0.4019,-0.6119,-0.493,-0.677,-0.04947,-0.01558,1.923,0.6792,2.411,0.3451,0.8794,0.5365,-0.2324,1.365,1.366,1.206,-0.0499,0.6472,-0.03535,0.7507,0.3708,-2.177,-0.7241,-0.6399,0.7127,-0.6908,0.6294,0.2391,-0.4224,0.3271,-1.508,-0.5504,1.211,0.7903,1.156,0.03053,0.3273,-0.3993,-0.5463,-0.2255,3.104,-0.5586,-0.7567,0.9735,1.667,0.09834,-0.2487,0.1991,1.511,1.519,0.6549,0.4527,-1.403,-1.11,0.7453,-0.5812,-1.534,0.002968,1.293,0.3523,-0.8455,0.3854,1.022,-0.4519,-2.662,0.4393,0.6411,2.198,1.475,1.46,-1.792,1.132,-0.8713,-1.863,-0.9052,-1.844,-0.4432,-0.6874,1.081,-0.08727,0.1708,-0.6959,-0.518,0.3455,-0.1922,0.7749,1.992,0.6572,0.8746,0.8947,0.7477,-1.291,-0.978,-0.4779,-0.2567,0.216,-0.8829,-1.665,-0.2463,1.647,0.3375,0.8522,0.8609,0.3641,0.3342,1.65,1.065,-0.4685,-0.5479,1.137,0.3566,-0.2158,-0.3787,-1.844,1.042,-0.06045,0.1424,0.5531,0.5725,-1.115,-0.5686,-0.1196,0.921,0.4022,-2.449,0.5892,0.001498,-0.5798,2.476,-0.9597,0.9117,-0.5212,1.638,-0.3057,0.3543,-1.205,-0.3509,0.8639,-0.547,-0.2219,2.277,-0.1604,-0.156,1.653,-0.2542,1.021,-0.2926,-0.9607,0.5775,0.3705,-0.712,-0.5594,-1.1,0.09135,-1.072,1.1,0.6615,0.2534,-0.1757,0.0357,-0.4831,-1.648,1.444,1.485,0.7318,-0.004524,0.9932,2.489,1.328,0.4638,0.6462,0.6259,-1.384,0.6016,1.261,-0.2041,0.82,-0.2122,-1.283,0.8762,0.7938,-0.7485,0.08108,-1.517,0.168,-1.493,1.018,-0.6016,1.57,2.279,0.2193,1.151,-1.107,0.1768,0.08331,1.887,0.8434,-0.5591,0.1062,-0.7653,-0.1668,-0.6109,0.6428,0.5855,0.6023,-0.1984,-0.7815,0.7557,-3.525,-0.5938,0.3517,0.8686),
                                        nr = 5, nc = 49), "kj")
  formdettimecov <- list(lambda0 ~ dettime, sigma ~ 1, D ~ 1)
  mod <- suppressWarnings(ScrModel$new(formdettimecov, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.7574,0.07993,2.739,5.653,0.2885,0.1779,0.1149,0.3412,-1.323,-0.2688,2.514,4.985,-0.192,0.4287,2.964,6.322),
                                                      nr = 4, nc = 4), check.attributes = FALSE)
})

test_that("Detection covariate smooths work", {
  formsm <- list(lambda0 ~ s(dettime, k = 3), sigma ~ 1, D ~ 1)
  mod <- suppressWarnings(ScrModel$new(formsm, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.7466,0.09975,0.08802,2.739,5.654,0.2859,0.8272,0.1941,0.1149,0.3412,-1.307,-1.521,-0.2924,2.514,4.985,-0.1863,1.721,0.4684,2.964,6.322),
                                                      nr = 5, nc = 4), check.attributes = FALSE)
  
})

test_that("Density surface models work", {
  formd <- list(lambda0  ~ 1, sigma ~ 1, D ~ x*y)
  mod <- suppressWarnings(ScrModel$new(formd, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.7255,2.73,5.311,-0.8221,-0.7311,1.067,0.2857,0.115,0.4402,0.7154,0.7186,1.223,-1.286,2.505,4.448,-2.224,-2.14,-1.33,-0.1654,2.956,6.174,0.58,0.6773,3.464),
                                                      nr = 6, nc = 4), check.attributes = FALSE)
  
  formsd <- list(lambda0  ~ 1, sigma ~ 1, D ~ s(x, k = 3) + s(y, k = 3))
  mod <- suppressWarnings(ScrModel$new(formsd, scrdat, start, print = FALSE))
  suppressWarnings(mod$fit())
  expect_equal(signif(mod$estimates()$par, 4), matrix(c(-0.6743,2.715,6.753,-5.744,-0.6992,-3.376,-0.797,0.2844,0.1133,0.9068,3.537,0.5342,3.872,0.6217,-1.232,2.493,4.976,-12.68,-1.746,-10.96,-2.015,-0.1168,2.937,8.531,1.188,0.3478,4.213,0.4214),
                                                      nr = 7, nc = 4), check.attributes = FALSE)
}
)


## Useful code for testers 
# ch <- scrdat$capthist()
# covariates(traps(ch)) <- data.frame(detcov = c(-0.3625,-0.4869,1.678,0.2306,-1.673,-1.118,0.4265,
#                                                -2.006,0.5957,-0.4685,0.823,-1.23,-0.2205,1.673,
#                                                1.341,-1.483,0.2875,-1.081,1.014,-0.3805,1.369,
#                                                -1.407,-0.4858,-0.1298,2.413,2.015,1.099,
#                                                -1.879,-0.1079,2.272,-1.377,0.5986,-0.6217,
#                                                0.2381,-0.6856,-1.245,-0.09166,0.03824,-0.7979,
#                                                -0.2653,0.6487,-0.09794,-0.3511,-0.8161,0.5126,
#                                                0.7343,-0.5492,-1.074,-0.2868))
# covariates(traps(ch)) <- data.frame(dettime = t(matrix(c(-1.125,-0.3552,-0.7712,0.06758,1.097,0.4569,-0.9512,0.1436,1.354,0.9711,-0.4308,-0.1576,-0.3884,-0.0782,-0.9711,0.6355,0.7569,1.606,0.0005856,0.5376,0.9206,-0.0548,0.9233,0.9337,-0.8067,0.05913,-0.6932,1.272,-0.4794,0.0427,0.6546,0.02216,-1.147,0.4019,-0.6119,-0.493,-0.677,-0.04947,-0.01558,1.923,0.6792,2.411,0.3451,0.8794,0.5365,-0.2324,1.365,1.366,1.206,-0.0499,0.6472,-0.03535,0.7507,0.3708,-2.177,-0.7241,-0.6399,0.7127,-0.6908,0.6294,0.2391,-0.4224,0.3271,-1.508,-0.5504,1.211,0.7903,1.156,0.03053,0.3273,-0.3993,-0.5463,-0.2255,3.104,-0.5586,-0.7567,0.9735,1.667,0.09834,-0.2487,0.1991,1.511,1.519,0.6549,0.4527,-1.403,-1.11,0.7453,-0.5812,-1.534,0.002968,1.293,0.3523,-0.8455,0.3854,1.022,-0.4519,-2.662,0.4393,0.6411,2.198,1.475,1.46,-1.792,1.132,-0.8713,-1.863,-0.9052,-1.844,-0.4432,-0.6874,1.081,-0.08727,0.1708,-0.6959,-0.518,0.3455,-0.1922,0.7749,1.992,0.6572,0.8746,0.8947,0.7477,-1.291,-0.978,-0.4779,-0.2567,0.216,-0.8829,-1.665,-0.2463,1.647,0.3375,0.8522,0.8609,0.3641,0.3342,1.65,1.065,-0.4685,-0.5479,1.137,0.3566,-0.2158,-0.3787,-1.844,1.042,-0.06045,0.1424,0.5531,0.5725,-1.115,-0.5686,-0.1196,0.921,0.4022,-2.449,0.5892,0.001498,-0.5798,2.476,-0.9597,0.9117,-0.5212,1.638,-0.3057,0.3543,-1.205,-0.3509,0.8639,-0.547,-0.2219,2.277,-0.1604,-0.156,1.653,-0.2542,1.021,-0.2926,-0.9607,0.5775,0.3705,-0.712,-0.5594,-1.1,0.09135,-1.072,1.1,0.6615,0.2534,-0.1757,0.0357,-0.4831,-1.648,1.444,1.485,0.7318,-0.004524,0.9932,2.489,1.328,0.4638,0.6462,0.6259,-1.384,0.6016,1.261,-0.2041,0.82,-0.2122,-1.283,0.8762,0.7938,-0.7485,0.08108,-1.517,0.168,-1.493,1.018,-0.6016,1.57,2.279,0.2193,1.151,-1.107,0.1768,0.08331,1.887,0.8434,-0.5591,0.1062,-0.7653,-0.1668,-0.6109,0.6428,0.5855,0.6023,-0.1984,-0.7815,0.7557,-3.525,-0.5938,0.3517,0.8686),
#                                                      nr = 5, nc = 49)))
# timevaryingcov(traps(ch)) <- list(dettime = 1:5)
# tmod <- secr.fit(ch, mask = scrdat$mesh(), detectfn = "HHN", model = formsm, start = start)
# coef(tmod)
# 
# paste(as.vector(signif(mod$estimates()$par, 4)), collapse=",")



