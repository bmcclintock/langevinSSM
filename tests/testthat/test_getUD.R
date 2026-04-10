# tests/testthat/test_getUD.R

# --- Helpers ---
get_static_raster <- function(name = "habitat") {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10))
  terra::values(r) <- (1:100) / 100
  names(r) <- name
  return(r)
}

get_dynamic_raster <- function() {
  r1 <- get_static_raster("time1")
  r2 <- get_static_raster("time2")
  r_dyn <- c(r1, r2)
  terra::time(r_dyn) <- as.POSIXct(c("2023-01-01", "2023-01-02"), tz = "UTC")
  return(r_dyn)
}

get_mock_fit_for_ud <- function() {
  fit <- list(
    estimates = list(
      natural = matrix(c(0.5, -0.2), nrow = 2, dimnames = list(c("beta_hab1", "beta_hab2"), "Estimate")),
      cov_natural = matrix(c(0.1, 0.01, 0.01, 0.1), nrow = 2, dimnames = list(c("beta_hab1", "beta_hab2"), c("beta_hab1", "beta_hab2")))
    ),
    signatures = list(
      data = NULL,
      covs = NULL
    )
  )
  class(fit) <- "fitLangevin"
  return(fit)
}

# --- Base functionality tests ---

test_that("getUD enforces matching covariate and beta lengths", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster())

  # Provide 2 covariates but only 1 beta
  expect_error(getUD(spatialCovs = covs, beta = c(0.5)),
               "length\\(spatialCovs\\) must equal length\\(beta\\)")
})

test_that("getUD handles multiple covariates successfully", {
  covs_single <- list(hab1 = get_static_raster())
  covs_multi <- list(hab1 = get_static_raster(), hab2 = get_static_raster())

  expect_s4_class(getUD(covs_single, beta = 0.5), "SpatRaster")
  expect_s4_class(getUD(covs_multi, beta = c(0.5, -0.2)), "SpatRaster")
})

test_that("getUD computes static log UD correctly", {
  covs <- list(hab = get_static_raster())
  ud <- getUD(spatialCovs = covs, beta = c(1), log = TRUE)

  expect_s4_class(ud, "SpatRaster")
  expect_equal(terra::nlyr(ud), 1)
  expect_equal(names(ud), "log_UD")
})

test_that("getUD normalizes static probabilities when log = FALSE", {
  covs <- list(hab = get_static_raster())
  beta <- c(0.5)

  # Calculate raw log-UD and normalized probability UD
  ud_log <- getUD(spatialCovs = covs, beta = beta, log = TRUE)
  ud_prob <- getUD(spatialCovs = covs, beta = beta, log = FALSE)

  # The probability UD should sum exactly to 1
  sum_prob <- terra::global(ud_prob, "sum", na.rm = TRUE)$sum
  expect_equal(sum_prob, 1, tolerance = 1e-6)

  # The log-UD should NOT sum to 1
  sum_log <- terra::global(ud_log, "sum", na.rm = TRUE)$sum
  expect_false(isTRUE(all.equal(sum_log, 1)))

  expect_equal(names(ud_prob), "UD")
})

test_that("getUD properly formats dynamic (time-varying) log UDs", {
  covs_dyn <- list(hab = get_dynamic_raster())
  beta <- c(0.5)

  ud_dyn <- getUD(spatialCovs = covs_dyn, beta = beta, log = TRUE)

  # The resulting UD should have the same number of layers as the dynamic covariate
  expect_equal(terra::nlyr(ud_dyn), 2)

  # The time attribute should be preserved and passed to the output UD
  ud_times <- terra::time(ud_dyn)
  expect_false(is.null(ud_times))
  expect_equal(length(ud_times), 2)
  expect_true(all(grepl("^log_UD$", names(ud_dyn))))
})

test_that("getUD dynamic probability UD layers individually sum to 1", {
  covs_dyn <- list(hab = get_dynamic_raster())
  ud_prob_dyn <- getUD(spatialCovs = covs_dyn, beta = c(0.5), log = FALSE)

  sums_ud <- terra::global(ud_prob_dyn, "sum", na.rm = TRUE)$sum

  # Ensure normalization applies across all time slices
  expect_equal(sums_ud[1], 1, tolerance = 1e-6)
  expect_equal(sums_ud[2], 1, tolerance = 1e-6)
})

# --- Simulation and Uncertainty Tests (nsims > 0) ---

# --- Simulation and Uncertainty Tests (nsims > 0) ---

test_that("getUD marginal posterior simulation returns correctly structured list", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster("hab2"))
  fit <- get_mock_fit_for_ud()

  # Set show_progress = FALSE to keep test output clean
  out <- getUD(spatialCovs = covs, fit = fit, nsims = 5, log = FALSE, show_progress = FALSE)

  expect_type(out, "list")
  expect_named(out, c("UD", "SE", "CV"))
  expect_s4_class(out$UD, "SpatRaster")
  expect_s4_class(out$SE, "SpatRaster")
  expect_s4_class(out$CV, "SpatRaster")

  expect_equal(names(out$UD), "UD")
  expect_equal(names(out$SE), "UD_SE")
  expect_equal(names(out$CV), "UD_CV")
})

test_that("getUD simulation handles dynamic (time-varying) covariates properly", {
  covs_dyn <- list(hab1 = get_dynamic_raster())
  fit_dyn <- list(
    estimates = list(
      natural = matrix(c(0.5), nrow = 1, dimnames = list(c("beta_hab1"), "Estimate")),
      cov_natural = matrix(c(0.1), nrow = 1, dimnames = list(c("beta_hab1"), c("beta_hab1")))
    ),
    signatures = list(data = NULL, covs = NULL)
  )
  class(fit_dyn) <- "fitLangevin"

  out <- getUD(spatialCovs = covs_dyn, fit = fit_dyn, nsims = 5, show_progress = FALSE)

  # Each output metric should match the number of layers in the dynamic covariate
  expect_equal(terra::nlyr(out$UD), 2)
  expect_equal(terra::nlyr(out$SE), 2)
  expect_equal(terra::nlyr(out$CV), 2)

  # Time attributes must be passed to the uncertainty metrics
  expect_equal(length(terra::time(out$SE)), 2)
  expect_equal(length(terra::time(out$CV)), 2)
})

test_that("getUD catches user errors related to missing fit requirements for nsims > 0", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster("hab2"))
  fit_bad <- get_mock_fit_for_ud()
  fit_bad$estimates$cov_natural <- NULL

  # Providing beta instead of fit when nsims > 0
  expect_error(
    getUD(spatialCovs = covs, beta = c(0.5, -0.2), nsims = 5, show_progress = FALSE),
    "argument \"fit\" is missing"
  )

  # Providing a fit object that lacks the covariance matrix
  expect_error(
    getUD(spatialCovs = covs, fit = fit_bad, nsims = 5, show_progress = FALSE),
    "Refit model to get cov_natural"
  )
})
