# tests/testthat/test_lambda_estimate.R
#
# Round 7: par$lambda (natural scale) for lambda = "estimate".
# Verifies that the barrier strength can be supplied on the natural scale through
# initialValues / checkPar / fitLangevin$par, is reported on the natural scale by
# getPar / fit$conditions, is invariant to scaleFactor, and is rejected whenever
# lambda is not "estimate".

# The 5-location track is too short for the (default) variogram initialisation, which now returns an
# error, so initialValues() is called with type = "neutral" (fitLangevin() calls supply sigma and gamma).
make_barrier_covs <- function() {
  r <- terra::rast(nrows = 20, ncols = 20, ext = c(-100, 200, -100, 200))
  terra::values(r) <- ifelse(terra::xFromCell(r, seq_len(terra::ncell(r))) >= 50, 1, 0)
  names(r) <- "coast_barrier"
  bar <- suppressMessages(prepBarrier(r))
  cov1 <- terra::rast(nrows = 20, ncols = 20, ext = c(-100, 200, -100, 200))
  terra::values(cov1) <- terra::yFromCell(cov1, seq_len(terra::ncell(cov1))) / 100
  names(cov1) <- "cov1"
  list(cov1 = cov1, coast_barrier = bar)
}

make_barrier_data <- function() {
  df <- data.frame(
    id = as.factor(rep("A", 5)),
    date = as.POSIXct("2024-01-01 12:00:00", tz = "UTC") + c(0, 1, 2, 4, 5) * 3600,
    dt = c(0, 1, 1, 2, 1),
    x = c(60, 60, 30, 60, 60), y = c(50, 55, 60, 65, 70),
    x.err = c(1, 1, 20, 1, 1), y.err = c(1, 1, 20, 1, 1),
    smaj = NA, smin = NA, eor = NA, lc = as.factor("G"))
  attr(df, "time.unit") <- "hours"
  class(df) <- c("dataLangevin", "data.frame")
  df
}

fit_estimate <- function(dat, covs, lambda0 = 0.5, ...) {
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA))
  suppressWarnings(suppressMessages(fitLangevin(
    data = dat, model = "underdamped", spatialCovs = covs, barrier = "coast_barrier",
    par = list(beta = 0.1, sigma = 5, gamma = 0.5, lambda = lambda0),
    map = m, lambda = "estimate", silent = TRUE, calcSE = FALSE,
    control = list(trace = 0, iter.max = 100, eval.max = 100), ...)))
}

test_that("initialValues returns a natural-scale par$lambda for lambda = 'estimate'", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()

  iv <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral", lambda = "estimate"))
  expect_true("lambda" %in% names(iv))
  expect_true(is.numeric(iv$lambda) && length(iv$lambda) == 1 && iv$lambda > 0)
  # default start is the tied stability limit at the median step
  expect_equal(iv$lambda,
               langevinSSM:::.stability_limit("underdamped", iv$sigma, iv$gamma, 1),
               tolerance = 1e-8)

  # no lambda by default
  iv0 <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral"))
  expect_null(iv0$lambda)

  # a user-supplied start is respected
  iv1 <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral",
                                        par = list(lambda = 3), lambda = "estimate"))
  expect_equal(iv1$lambda, 3)
})

test_that("initialValues rejects par$lambda unless lambda = 'estimate'", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  expect_error(
    suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, par = list(lambda = 5))),
    "Cannot specify par\\$lambda")
})

test_that("checkPar converts par$lambda to a free log_lambda (working scale)", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  iv <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral", lambda = "estimate"))
  iv$lambda <- 2.5

  cp <- langevinSSM:::checkPar(iv, "underdamped", spatialCovs = covs, lambda_mode = 2L)
  expect_equal(cp$par$log_lambda, log(2.5), tolerance = 1e-12)
  expect_null(cp$par$lambda)
  expect_true(is.null(cp$map$log_lambda))          # free by default

  # map$lambda is moved to map$log_lambda
  cp2 <- langevinSSM:::checkPar(iv, "underdamped", spatialCovs = covs,
                                map = list(lambda = factor(NA)), lambda_mode = 2L)
  expect_true(is.na(cp2$map$log_lambda))
})

test_that("checkPar rejects 'lambda' in par/map when lambda_mode != 2", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  iv <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral", lambda = "estimate"))
  expect_error(langevinSSM:::checkPar(iv, "underdamped", spatialCovs = covs, lambda_mode = 0L),
               "can only be specified")
  expect_error(langevinSSM:::checkPar(iv, "underdamped", spatialCovs = covs, lambda_mode = 1L),
               "can only be specified")
})

test_that("checkPar requires par$lambda when lambda = 'estimate'", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  iv <- suppressMessages(initialValues(dat, model = "underdamped", spatialCovs = covs, type = "neutral", lambda = "estimate"))
  iv$lambda <- NULL
  expect_error(langevinSSM:::checkPar(iv, "underdamped", spatialCovs = covs, lambda_mode = 2L),
               "par\\$lambda is missing")
})

test_that("fitLangevin(lambda = 'estimate') accepts par$lambda and reports it on the natural scale", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  fit <- fit_estimate(dat, covs, lambda0 = 0.5)

  expect_identical(fit$conditions$lambda_mode, "estimate")
  expect_true(is.finite(fit$conditions$lambda) && fit$conditions$lambda > 0)
  expect_equal(getPar(fit)$lambda, fit$conditions$lambda, tolerance = 1e-8)
  expect_true("lambda" %in% rownames(fit$estimates$natural))
  expect_equal(fit$estimates$natural["lambda", "Estimate"], fit$conditions$lambda, tolerance = 1e-8)
})

test_that("the reported natural-scale lambda is invariant to scaleFactor", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA))

  lr <- vapply(c(1, 2), function(sf) {
    fit <- suppressWarnings(suppressMessages(fitLangevin(
      data = dat, model = "underdamped", spatialCovs = covs, barrier = "coast_barrier",
      par = list(beta = 0.1, sigma = 5, gamma = 0.5), map = m, lambda = "tied",
      scaleFactor = sf, silent = TRUE, calcSE = FALSE)))
    fit$conditions$lambda_report
  }, numeric(1))

  expect_equal(lr[1], lr[2], tolerance = 1e-8)
})

test_that("a prior on log_lambda (natural scale) is applied when lambda = 'estimate'", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  pr <- data.frame(mean = 0, sd = 0.5); rownames(pr) <- "log_lambda"

  # prior mean 0 => lambda = 1, far from the weakly-identified optimum without a prior
  fit0 <- fit_estimate(dat, covs, lambda0 = 100)
  fit <- fit_estimate(dat, covs, lambda0 = 100, prior = pr)
  l0 <- getPar(fit0)$lambda; l1 <- getPar(fit)$lambda
  expect_true(abs(log(l1)) < abs(log(l0)))          # pulled towards the prior mean
  expect_true(abs(log(l1)) < 0.25)
  # a tight prior pins lambda at its mean
  pr$sd <- 0.01
  fit2 <- fit_estimate(dat, covs, lambda0 = 100, prior = pr)
  expect_equal(getPar(fit2)$lambda, 1, tolerance = 1e-2)
})

test_that("a log_lambda prior is rejected unless lambda = 'estimate'", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA))
  pr <- data.frame(mean = 0, sd = 0.5); rownames(pr) <- "log_lambda"
  expect_error(suppressMessages(fitLangevin(
    data = dat, model = "underdamped", spatialCovs = covs, barrier = "coast_barrier",
    par = list(beta = 0.1, sigma = 5, gamma = 0.5), map = m, lambda = "tied",
    prior = pr, silent = TRUE, calcSE = FALSE)), "log_lambda")
})

test_that("tied and fixed lambda reject par$lambda", {
  covs <- make_barrier_covs(); dat <- make_barrier_data()
  m <- list(beta = factor(NA), sigma = factor(NA), gamma = factor(NA))

  expect_error(suppressMessages(fitLangevin(
    data = dat, model = "underdamped", spatialCovs = covs, barrier = "coast_barrier",
    par = list(beta = 0.1, sigma = 5, gamma = 0.5, lambda = 100), map = m,
    lambda = "tied", silent = TRUE, calcSE = FALSE)), "par\\$lambda")

  expect_error(suppressMessages(fitLangevin(
    data = dat, model = "underdamped", spatialCovs = covs, barrier = "coast_barrier",
    par = list(beta = 0.1, sigma = 5, gamma = 0.5, lambda = 100), map = m,
    lambda = 100, silent = TRUE, calcSE = FALSE)), "par\\$lambda")
})
