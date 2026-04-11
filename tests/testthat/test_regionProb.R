# tests/testthat/test_regionProb.R

# --- Helpers ---
get_mock_covs_and_mask <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10))
  terra::values(r) <- runif(100)
  names(r) <- "habitat"

  # Create a mask covering half the raster
  mask <- r
  terra::values(mask) <- ifelse(1:100 > 50, 1, 0)

  return(list(covs = list(habitat = r), mask = mask))
}

get_mock_fit_for_rp <- function() {
  fit <- list(
    estimates = list(
      natural = data.frame(Estimate = c(0.5, 2.0), row.names = c("beta_habitat", "sigma"))
    ),
    covariance = list(
      natural = matrix(c(0.1, 0, 0, 0.05), nrow = 2, dimnames = list(c("beta_habitat", "sigma"), c("beta_habitat", "sigma")))
    ),
    signatures = list(data = NULL, covs = NULL)
  )
  class(fit) <- "fitLangevin"
  return(fit)
}

test_that("regionProb catches invalid inputs", {
  env <- get_mock_covs_and_mask()
  fit <- get_mock_fit_for_rp()

  # 1. Invalid Fit
  expect_error(regionProb(fit = list(), spatialCovs = env$covs, mask = env$mask),
               "must be a fitLangevin object")

  # 2. Invalid Spatial Covariates
  expect_error(regionProb(fit = fit, spatialCovs = env$mask, mask = env$mask),
               "must be a list of SpatRaster")

  # 3. Geometry mismatch
  bad_mask <- terra::rast(nrows = 5, ncols = 5, ext = c(0, 5, 0, 5))
  expect_error(regionProb(fit = fit, spatialCovs = env$covs, mask = bad_mask),
               "same projection \\(CRS\\), extent, and resolution")

  # 4. Invalid nSims
  expect_error(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = -10),
               "single non-negative integer")
})

test_that("regionProb correctly calculates Delta method and Simulation bounds", {
  env <- get_mock_covs_and_mask()
  fit <- get_mock_fit_for_rp()

  out <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 10, show_progress = FALSE))

  # Check Structure
  expect_type(out, "list")
  expect_named(out, c("Point_Estimate", "SE_delta", "CI_delta_95", "SE_sim", "CI_sim_95", "simulated_draws"))

  # Point Estimate Bounds
  expect_true(out$Point_Estimate > 0 && out$Point_Estimate < 1)

  # Delta Bounds Valid
  expect_true(out$CI_delta_95[1] >= 0)
  expect_true(out$CI_delta_95[2] <= 1)
  expect_true(out$CI_delta_95[1] <= out$CI_delta_95[2])

  # Sim Bounds Valid
  expect_length(out$simulated_draws, 10)
  expect_true(out$CI_sim_95[1] >= 0)
  expect_true(out$CI_sim_95[2] <= 1)
})

test_that("regionProb appropriately caps probability and SE when hitting limits", {
  # If the mask is entirely 1s, the probability should perfectly equal 1,
  # and the SE should be forced to 0
  env <- get_mock_covs_and_mask()
  fit <- get_mock_fit_for_rp()

  # Force mask to all 1s
  terra::values(env$mask) <- 1

  out <- regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 0, show_progress = FALSE)

  # Sum of all probabilities must be 1, so SE must be 0
  expect_equal(out$Point_Estimate, 1)
  expect_equal(out$SE_delta, 0)
  expect_equal(out$CI_delta_95, c(1, 1))
})
