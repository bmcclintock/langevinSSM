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

test_that("getUD enforces matching covariate and beta lengths", {
  covs <- list(hab1 = get_static_raster(), hab2 = get_static_raster())

  # Provide 2 covariates but only 1 beta
  expect_error(getUD(spatialCovs = covs, beta = c(0.5)),
               "length\\(spatialCovs\\) must equal length\\(beta\\)")
})

test_that("getUD normalizes probabilities when log = FALSE", {
  covs <- list(hab = get_static_raster())
  beta <- c(0.5)

  # Calculate raw log-UD
  ud_log <- getUD(spatialCovs = covs, beta = beta, log = TRUE)

  # Calculate normalized probability UD
  ud_prob <- getUD(spatialCovs = covs, beta = beta, log = FALSE)

  # 1. The probability UD should sum exactly to 1 (accounting for tiny floating point differences)
  sum_prob <- terra::global(ud_prob, "sum", na.rm = TRUE)$sum
  expect_equal(sum_prob, 1, tolerance = 1e-6)

  # 2. The log-UD should NOT sum to 1
  sum_log <- terra::global(ud_log, "sum", na.rm = TRUE)$sum
  expect_false(isTRUE(all.equal(sum_log, 1)))
})

test_that("getUD handles multiple covariates successfully", {
  # This explicitly tests the fix we made earlier to avoid the 2:length(x) bug
  covs_single <- list(hab1 = get_static_raster())
  covs_multi <- list(hab1 = get_static_raster(), hab2 = get_static_raster())

  expect_s4_class(getUD(covs_single, beta = 0.5), "SpatRaster")
  expect_s4_class(getUD(covs_multi, beta = c(0.5, -0.2)), "SpatRaster")
})

test_that("getUD properly formats dynamic (time-varying) rasters", {
  covs_dyn <- list(hab = get_dynamic_raster())
  beta <- c(0.5)

  ud_dyn <- getUD(spatialCovs = covs_dyn, beta = beta)

  # The resulting UD should have the same number of layers as the dynamic covariate
  expect_equal(terra::nlyr(ud_dyn), 2)

  # The time attribute should be preserved and passed to the output UD
  ud_times <- terra::time(ud_dyn)
  expect_false(is.null(ud_times))
  expect_equal(length(ud_times), 2)
})
