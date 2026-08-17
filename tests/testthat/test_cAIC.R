# tests/testthat/test_cAIC.R

# --- Setup Shared Test Data ---
# We use a very small track to ensure rapid execution during CRAN checks
set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
testPar <- list(beta = c(-1, 1, 0, 0), sigma = 2, gamma = 0.5)
testMeas <- list(smaj.sd = 1.5, smin.sd = 0.75, eor.lim = c(0, 180))

# Simulate a micro-track using the package's built-in exampleCovs
testDat <- suppressWarnings(suppressMessages(simLangevin(
  par = testPar,
  spatialCovs = exampleCovs,
  nbAnimals = 1,
  obsPerAnimal = 30,
  measurementError = testMeas
)))

test_that("cAIC strictly enforces input object classes", {
  bad_fit <- list(par = c(1, 2))
  class(bad_fit) <- "lm"

  expect_error(
    cAIC(bad_fit, testDat, exampleCovs),
    regexp = "'fit' must be a 'fitLangevin' object"
  )
})

test_that("cAIC catches missing joint precision matrices (Trap #1)", {
  # Fit a model explicitly turning OFF the joint precision matrix.
  # suppressWarnings silences the expected convergence warnings caused by iter.max = 0.
  fit_no_prec <- suppressWarnings(suppressMessages(fitLangevin(
    data = testDat,
    spatialCovs = exampleCovs,
    model = "overdamped",
    getJointPrecision = FALSE,
    silent = TRUE,
    control = list(eval.max = 0, iter.max = 0)
  )))

  expect_error(
    cAIC(fit_no_prec, testDat, exampleCovs),
    regexp = "The model must be fitted with 'getJointPrecision = TRUE'"
  )
})

test_that("cAIC triggers safeguards upon signature mismatches (Trap #2)", {
  # Perform a valid fit to ensure jointPrecision extraction passes without NaNs
  fit_base <- suppressWarnings(suppressMessages(fitLangevin(
    data = testDat,
    spatialCovs = exampleCovs,
    model = "overdamped",
    getJointPrecision = TRUE,
    silent = TRUE
  )))

  # Trap 2A: The user tampers with the data frame prior to cAIC evaluation
  tampered_dat <- testDat
  tampered_dat$x[1] <- tampered_dat$x[1] + 500

  expect_error(
    cAIC(fit_base, tampered_dat, exampleCovs),
    regexp = "Safeguard triggered: the 'data' provided does not match"
  )

  # Trap 2B: The user tampers with the spatial covariates
  tampered_covs <- exampleCovs
  terra::values(tampered_covs[[1]]) <- terra::values(tampered_covs[[1]]) * 2

  expect_error(
    cAIC(fit_base, testDat, tampered_covs),
    regexp = "Safeguard triggered: the 'spatialCovs' provided do not match"
  )
})

test_that("cAIC correctly computes Method 2 for an Overdamped model", {
  fit_od <- suppressWarnings(suppressMessages(fitLangevin(
    data = testDat,
    spatialCovs = exampleCovs,
    model = "overdamped",
    getJointPrecision = TRUE,
    silent = TRUE
  )))

  res <- suppressWarnings(suppressMessages(cAIC(fit_od, testDat, exampleCovs)))

  expect_s3_class(res, "caicLangevin")
  expect_type(res$cAIC, "double")
  expect_type(res$trace_penalty, "double")

  # Check random effect dimensionality: 30 obs * 2 coordinates (mu.x, mu.y) = 60
  expect_equal(res$q, 60)

  # Ensure the Effective Degrees of Freedom incorporates the trace penalty
  expect_equal(res$EDF, res$q - res$trace_penalty)
})

test_that("cAIC correctly computes Method 2 for an Underdamped model", {
  fit_ud <- suppressWarnings(suppressMessages(fitLangevin(
    data = testDat,
    spatialCovs = exampleCovs,
    model = "underdamped",
    getJointPrecision = TRUE,
    silent = TRUE
  )))

  res <- suppressWarnings(suppressMessages(cAIC(fit_ud, testDat, exampleCovs)))

  # Check random effect dimensionality: 30 obs * 2 coords (mu) + 30 obs * 2 velocities (vel) = 120
  expect_equal(res$q, 120)

  # Conditional NLL must be finite and computable
  expect_true(is.finite(res$NLL_cond))
  expect_true(res$trace_penalty > 0)
})

test_that("Print method successfully executes", {
  # A proper fit ensures the Hessian is positive-definite, preventing sdreport crashes
  fit_od <- suppressWarnings(suppressMessages(fitLangevin(
    data = testDat,
    spatialCovs = exampleCovs,
    model = "overdamped",
    getJointPrecision = TRUE,
    silent = TRUE
  )))

  res <- suppressWarnings(suppressMessages(cAIC(fit_od, testDat, exampleCovs)))

  # Test that the S3 print method properly registers and doesn't throw formatting errors
  expect_output(print(res), "Conditional Akaike Information Criterion")
  expect_output(print(res), "Trace Penalty:")
})
