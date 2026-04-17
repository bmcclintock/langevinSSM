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
      natural = matrix(c(0.5, -0.2), nrow = 2, dimnames = list(c("beta_hab1", "beta_hab2"), "Estimate"))
    ),
    covariance = list(
      natural = matrix(c(0.1, 0.01, 0.01, 0.1), nrow = 2, dimnames = list(c("beta_hab1", "beta_hab2"), c("beta_hab1", "beta_hab2")))
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
  expect_error(getUD(spatialCovs = covs, beta = c(0.5), plot = FALSE), "length\\(spatialCovs\\) must equal length\\(beta\\)")
})

test_that("getUD handles multiple covariates successfully", {
  covs_single <- list(hab1 = get_static_raster())
  covs_multi <- list(hab1 = get_static_raster(), hab2 = get_static_raster())

  out_single <- getUD(covs_single, beta = 0.5, plot = FALSE)
  out_multi <- getUD(covs_multi, beta = c(0.5, -0.2), plot = FALSE)

  expect_s4_class(out_single, "SpatRaster")
  expect_s4_class(out_multi, "SpatRaster")
})

test_that("getUD computes static log UD correctly", {
  covs <- list(hab = get_static_raster())
  ud <- getUD(spatialCovs = covs, beta = c(1), log = TRUE, plot = FALSE)

  expect_s4_class(ud, "SpatRaster")
  expect_equal(terra::nlyr(ud), 1)
  expect_equal(names(ud), "log_UD")
})

test_that("getUD normalizes static probabilities when log = FALSE", {
  covs <- list(hab = get_static_raster())
  beta <- c(0.5)

  ud_log <- getUD(spatialCovs = covs, beta = beta, log = TRUE, plot = FALSE)
  ud_prob <- getUD(spatialCovs = covs, beta = beta, log = FALSE, plot = FALSE)

  sum_prob <- terra::global(ud_prob[["UD"]], "sum", na.rm = TRUE)$sum
  expect_equal(sum_prob, 1, tolerance = 1e-6)

  sum_log <- terra::global(ud_log[["log_UD"]], "sum", na.rm = TRUE)$sum
  expect_false(isTRUE(all.equal(sum_log, 1)))

  expect_equal(names(ud_prob), "UD")
})

test_that("getUD properly formats dynamic (time-varying) log UDs", {
  covs_dyn <- list(hab = get_dynamic_raster())
  ud_dyn <- getUD(spatialCovs = covs_dyn, beta = c(0.5), log = TRUE, plot = FALSE)

  expect_equal(terra::nlyr(ud_dyn), 2)
  ud_times <- terra::time(ud_dyn)
  expect_false(is.null(ud_times))
  expect_equal(length(ud_times), 2)
  expect_true(all(grepl("^log_UD$", names(ud_dyn))))
})

test_that("getUD dynamic probability UD layers individually sum to 1", {
  covs_dyn <- list(hab = get_dynamic_raster())
  ud_prob_dyn <- getUD(spatialCovs = covs_dyn, beta = c(0.5), log = FALSE, plot = FALSE)

  sums_ud <- terra::global(ud_prob_dyn, "sum", na.rm = TRUE)$sum
  expect_equal(sums_ud[1], 1, tolerance = 1e-6)
  expect_equal(sums_ud[2], 1, tolerance = 1e-6)
})

# --- Simulation and Uncertainty Tests (nSims > 0) ---

test_that("getUD marginal Monte Carlo simulation returns correctly structured raster stack", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster("hab2"))
  fit <- get_mock_fit_for_ud()

  out <- suppressMessages(getUD(spatialCovs = covs, fit = fit, nSims = 5, log = FALSE, show_progress = FALSE, plot = FALSE))

  expect_s4_class(out, "SpatRaster")
  expect_equal(terra::nlyr(out), 5)
  expected_names <- c("UD", "UD_SE_delta", "UD_CV_delta", "UD_SE_sim", "UD_CV_sim")
  expect_true(all(expected_names %in% names(out)))
})

test_that("getUD simulation handles dynamic (time-varying) covariates properly", {
  covs_dyn <- list(hab1 = get_dynamic_raster())
  fit_dyn <- list(
    estimates = list(natural = matrix(c(0.5), nrow = 1, dimnames = list(c("beta_hab1"), "Estimate"))),
    covariance = list(natural = matrix(c(0.1), nrow = 1, dimnames = list(c("beta_hab1"), c("beta_hab1")))),
    signatures = list(data = NULL, covs = NULL)
  )
  class(fit_dyn) <- "fitLangevin"

  out <- suppressMessages(getUD(spatialCovs = covs_dyn, fit = fit_dyn, nSims = 5, show_progress = FALSE, plot = FALSE))

  expect_equal(terra::nlyr(out), 10)
  expect_equal(length(terra::time(out)), 10)
  expect_equal(sum(names(out) == "UD_SE_delta"), 2)
})

test_that("getUD catches user errors related to missing fit requirements for nSims > 0", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster("hab2"))
  fit_bad <- get_mock_fit_for_ud()
  fit_bad$covariance$natural <- NULL

  expect_error(getUD(spatialCovs = covs, beta = c(0.5, -0.2), nSims = 5, show_progress = FALSE, plot = FALSE), "fit\\$covariance\\$natural not found")
  expect_error(getUD(spatialCovs = covs, fit = fit_bad, nSims = 5, show_progress = FALSE, plot = FALSE), "Refit model to get covariance matrix.")
})
