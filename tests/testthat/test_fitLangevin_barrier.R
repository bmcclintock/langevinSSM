# tests/testthat/test_fitLangevin_barrier.R

# ==============================================================================
# 1. SETUP MOCK DATA
# ==============================================================================

# Create a small 10x10 dummy barrier raster
r_ext <- ext(0, 10, 0, 10)
mock_rast <- rast(nrows = 10, ncols = 10, ext = r_ext)
terra::values(mock_rast) <- rep(c(0, 1), each = 50) # Bottom half water (0), Top half land (1)
names(mock_rast) <- "coast_mask"
mock_covs <- list(coast_mask = suppressMessages(prepBarrier(mock_rast)))

# Create a valid (but tiny) mock dataLangevin object
mock_df <- data.frame(
  id = factor(rep(1, 6)),
  x = c(5, 5.1, 5.2, 5.1, 5.0, 4.9),
  y = c(6, 6.1, 6.0, 5.9, 6.0, 6.1),
  date = as.POSIXct("2024-01-01 00:00:00") + seq(0, 5 * 3600, by = 3600),
  dt = c(0, 3600, 3600, 3600, 3600, 3600)
)
class(mock_df) <- c("dataLangevin", "data.frame")
attr(mock_df, "time.unit") <- "secs"

# ==============================================================================
# 2. TEST INPUT VALIDATION & USER ERRORS
# ==============================================================================

test_that("fitLangevin_barrier fails fast on invalid data classes", {
  bad_df <- mock_df
  class(bad_df) <- "data.frame"

  expect_error(
    fitLangevin_barrier(bad_df, mock_covs),
    "'data' must be a dataLangevin object"
  )
})

test_that("fitLangevin_barrier fails fast on insufficient data for empirical calculations", {
  tiny_df <- mock_df[1:3, ]
  class(tiny_df) <- c("dataLangevin", "data.frame")

  expect_error(
    suppressMessages(fitLangevin_barrier(tiny_df, mock_covs, lambda_max = NULL)),
    "Insufficient valid steps"
  )
})

test_that("fitLangevin_barrier handles spatial overlap failure gracefully", {
  oob_df <- mock_df
  oob_df$x <- oob_df$x + 500
  class(oob_df) <- c("dataLangevin", "data.frame")

  expect_error(
    suppressMessages(fitLangevin_barrier(oob_df, mock_covs, lambda_max = NULL)),
    "No valid spatial observations"
  )
})

# ==============================================================================
# 3. TEST EDGE CASES & LOGIC OVERRIDES
# ==============================================================================

test_that("fitLangevin_barrier strictly enforces minimum data size even if lambda_max is provided", {
  tiny_df <- mock_df[1:3, ]
  class(tiny_df) <- c("dataLangevin", "data.frame")

  # Because empirical_sigma is required for lambda_min, the function
  # must fail fast regardless of whether lambda_max is supplied.
  expect_error(
    suppressMessages(fitLangevin_barrier(tiny_df, mock_covs, lambda_max = 50, n_coarse = 2, n_sims = 1)),
    regexp = "Insufficient valid steps"
  )
})

test_that("fitLangevin_barrier handles total optimization failure gracefully", {
  # Use the valid mock_covs, but set lambda_max absurdly high to break numerical integration
  expect_error(
    suppressMessages(fitLangevin_barrier(mock_df, mock_covs,
                                         lambda_max = 1e6, n_coarse = 2, n_fine = 1, n_sims = 1)),
    "All models in the coarse grid failed"
  )
})

test_that("Initialization grids calculate safely without throwing internal NA/NaN warnings", {
  expect_warning(
    tryCatch(
      suppressMessages(fitLangevin_barrier(mock_df, mock_covs,
                                           n_coarse = 2, n_fine = 1, n_sims = 1)),
      error = function(e) NA
    ),
    regexp = NA
  )
})

test_that("fitLangevin_barrier fails fast if no barrier attribute is found", {
  # Strip the 'barLangevin' attribute from the mock covariates
  untagged_covs <- mock_covs
  attr(untagged_covs$coast_mask, "barLangevin") <- NULL

  expect_error(
    fitLangevin_barrier(mock_df, untagged_covs),
    "No barrier found in 'spatialCovs'. Did you run prepBarrier"
  )
})

test_that("fitLangevin_barrier fails fast if multiple barriers are tagged", {
  # Duplicate the tagged barrier raster
  multi_covs <- mock_covs
  multi_covs$second_mask <- mock_covs$coast_mask

  expect_error(
    fitLangevin_barrier(mock_df, multi_covs),
    "Multiple rasters in 'spatialCovs' are tagged as barriers"
  )
})
