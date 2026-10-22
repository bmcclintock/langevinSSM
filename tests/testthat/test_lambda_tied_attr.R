library(testthat)
library(langevinSSM)

# simulated coast: land for x < -600
tied_covs <- function() {
  r <- terra::rast(nrows = 40, ncols = 40, ext = c(-1000, 1000, -1000, 1000))
  terra::values(r) <- ifelse(terra::xFromCell(r, seq_len(terra::ncell(r))) >= -600, 1, 0)
  names(r) <- "coast_barrier"
  bar <- suppressMessages(prepBarrier(r))
  hab <- terra::rast(nrows = 40, ncols = 40, ext = c(-1000, 1000, -1000, 1000))
  terra::values(hab) <- terra::yFromCell(hab, seq_len(terra::ncell(hab))) / 1000
  names(hab) <- "habitat"
  list(habitat = hab, coast_barrier = bar)
}
tied_par <- list(beta = 1, sigma = 5, gamma = 0.5, psi = 1)
tied_me <- list(smaj.sd = 1.5, smin.sd = 0.75)

test_that("tied simLangevin stores the true lambda in attr 'lambda_tied' and keeps attr 'lambda' = 'tied'", {
  covs <- tied_covs()
  for (model in c("underdamped", "overdamped")) {
    p <- if (model == "underdamped") tied_par else tied_par[c("beta", "sigma", "psi")]
    g <- if (model == "underdamped") 0.5 else 1
    set.seed(1)
    sim <- suppressMessages(simLangevin(model = model, par = p, spatialCovs = covs, barrier = "coast_barrier",
                                        obsPerAnimal = 50, timeStep = 0.1, initialPosition = c(-500, 0), measurementError = tied_me))
    expect_identical(attr(sim, "lambda"), "tied")
    expect_equal(attr(sim, "lambda_tied"), langevinSSM:::.stability_limit(model, 5, g, 0.1))
    if (model == "overdamped") expect_equal(attr(sim, "lambda_tied"), 2 / (25 * 0.1))
  }
})

test_that("'lambda_tied' refers to the observed steps after subSample", {
  covs <- tied_covs()
  set.seed(2)
  sim <- suppressMessages(simLangevin(model = "underdamped", par = tied_par, spatialCovs = covs, barrier = "coast_barrier",
                                      obsPerAnimal = 400, timeStep = 0.05, initialPosition = c(-500, 0), measurementError = tied_me,
                                      subSample = list(samplingRate = 10, propMissing = 0)))
  expect_identical(attr(sim, "lambda"), "tied")
  med <- median(sim$dt[sim$dt > 1e-6])   # subSampleData samples rows at random, so the steps vary
  expect_gt(med, 0.05)
  expect_equal(attr(sim, "lambda_tied"), langevinSSM:::.stability_limit("underdamped", 5, 0.5, med))
  # smaller than the value used at the simulation step
  expect_lt(attr(sim, "lambda_tied"), langevinSSM:::.stability_limit("underdamped", 5, 0.5, 0.05))
})

test_that("numeric lambda and no barrier give no 'lambda_tied' attribute", {
  covs <- tied_covs()
  set.seed(3)
  sim_num <- suppressMessages(simLangevin(model = "underdamped", par = tied_par, spatialCovs = covs, barrier = "coast_barrier",
                                          lambda = 2, obsPerAnimal = 20, initialPosition = c(-500, 0), measurementError = tied_me))
  expect_equal(attr(sim_num, "lambda"), 2)
  expect_null(attr(sim_num, "lambda_tied"))
  sim_nb <- suppressMessages(simLangevin(model = "underdamped", par = tied_par, spatialCovs = covs["habitat"],
                                         obsPerAnimal = 20, initialPosition = c(-500, 0), measurementError = tied_me))
  expect_null(attr(sim_nb, "lambda"))
  expect_null(attr(sim_nb, "lambda_tied"))
})

test_that("'lambda_tied' equals the lambda_report of a fit at the true sigma and gamma; downstream use is unchanged", {
  covs <- tied_covs()
  set.seed(4)
  sim <- suppressMessages(simLangevin(model = "underdamped", par = tied_par, spatialCovs = covs, barrier = "coast_barrier",
                                      obsPerAnimal = 300, timeStep = 0.05, initialPosition = c(-500, 0),
                                      measurementError = tied_me,
                                      subSample = list(samplingRate = 10, propMissing = 0)))
  # fit with all parameters fixed at the truth: lambda_report is then the true tied lambda
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA), psi = factor(NA))
  fit <- suppressWarnings(suppressMessages(fitLangevin(data = sim, model = "underdamped", spatialCovs = covs,
                          barrier = "coast_barrier", lambda = NULL, par = tied_par, map = m, silent = TRUE, calcSE = FALSE)))
  expect_identical(fit$conditions$lambda_mode, "tied")   # auto-detected from attr(, "lambda")
  expect_equal(attr(sim, "lambda_tied"), fit$conditions$lambda_report, tolerance = 1e-8)

  # getUD with the stored attribute is still a hard wall
  ud <- getUD(covs, beta = 1, barrier = "coast_barrier", lambda = attr(sim, "lambda"), log = FALSE, plot = FALSE)
  land <- terra::values(covs$coast_barrier)[, 1] <= 0
  expect_equal(sum(terra::values(ud)[land, 1], na.rm = TRUE), 0)

  # simLangevin from the tied fit stores the fit's lambda (fit$conditions$lambda_report)
  for (cnd in c(TRUE, FALSE)) {
    set.seed(5)
    s2 <- suppressWarnings(suppressMessages(simLangevin(fit, data = sim, spatialCovs = covs, timeStep = 0.05, conditional = cnd)))
    expect_identical(attr(s2, "lambda"), "tied")
    expect_identical(attr(s2, "lambda_tied"), fit$conditions$lambda_report)
  }
})

test_that("a fixed-lambda fit gives no 'lambda_tied' attribute from simLangevin", {
  covs <- tied_covs()
  set.seed(6)
  sim <- suppressMessages(simLangevin(model = "underdamped", par = tied_par, spatialCovs = covs, barrier = "coast_barrier",
                                      lambda = 0.5, obsPerAnimal = 100, timeStep = 0.5, initialPosition = c(-500, 0),
                                      measurementError = tied_me))
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA), psi = factor(NA))
  fit <- suppressWarnings(suppressMessages(fitLangevin(data = sim, model = "underdamped", spatialCovs = covs,
                          barrier = "coast_barrier", lambda = 0.5, par = tied_par, map = m, silent = TRUE, calcSE = FALSE)))
  s2 <- suppressWarnings(suppressMessages(simLangevin(fit, data = sim, spatialCovs = covs, timeStep = 0.5, conditional = FALSE)))
  expect_equal(attr(s2, "lambda"), 0.5)
  expect_null(attr(s2, "lambda_tied"))
})
