# tests/testthat/test_suggestLambda.R

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

test_that("suggestLambda fails fast on invalid data classes", {
  bad_df <- mock_df
  class(bad_df) <- "data.frame"

  expect_error(
    suggestLambda(bad_df, mock_covs, barrier = "coast_mask", timeStep = "1 hour"),
    "'data' must be a dataLangevin object"
  )
})

test_that("suggestLambda fails fast on insufficient data for empirical calculations", {
  tiny_df <- mock_df[1:3, ]
  class(tiny_df) <- c("dataLangevin", "data.frame")
  attr(tiny_df, "time.unit") <- "secs"

  expect_error(
    suppressMessages(suggestLambda(tiny_df, mock_covs, barrier = "coast_mask")),
    "At least one habitat covariate must"
  )
})

test_that("suggestLambda handles spatial overlap failure gracefully", {
  oob_df <- mock_df
  oob_df$x <- oob_df$x + 500
  class(oob_df) <- c("dataLangevin", "data.frame")
  attr(oob_df, "time.unit") <- "secs"

  expect_error(
    suppressMessages(suggestLambda(oob_df, spatialCovs=list(mock=mock_covs$coast_mask,coast_mask=mock_covs$coast_mask), barrier = "coast_mask")),
    "The tracking data do not overlap"
  )
})

# ==============================================================================
# 3. TEST EDGE CASES & LOGIC OVERRIDES
# ==============================================================================

test_that("suggestLambda fails fast if no barrier attribute is found", {
  # Strip the 'barLangevin' attribute from the mock covariates
  untagged_covs <- mock_covs
  attr(untagged_covs$coast_mask, "barLangevin") <- NULL

  expect_error(
    suppressMessages(suggestLambda(mock_df, untagged_covs, barrier = "coast_mask")),
    "barrier is not a 'barLangevin' object created by prepBarrier"
  )
})

test_that("suggestLambda fails fast if barrier name not found in spatialCovs", {
  # Pass a string for a raster that isn't in the list
  expect_error(
    suppressMessages(suggestLambda(mock_df, mock_covs, barrier = "missing_mask")),
    "Barrier raster 'missing_mask' not found in spatialCovs"
  )
})
