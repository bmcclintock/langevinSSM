# tests/testthat/test_rasterOverlap.R

# --- Helper functions for mock data ---
get_mock_raster <- function(vals = 1, extent = c(0, 10, 0, 10), res = 1) {
  r <- terra::rast(ext = extent, res = res)
  terra::values(r) <- rep(vals, length.out = terra::ncell(r))
  return(r)
}

test_that("rasterOverlap catches non-SpatRaster inputs", {
  r1 <- get_mock_raster()
  bad_input <- matrix(1, nrow = 10, ncol = 10)

  expect_error(rasterOverlap(r1, bad_input), "Both inputs must be terra::SpatRaster objects")
  expect_error(rasterOverlap(bad_input, r1), "Both inputs must be terra::SpatRaster objects")
})

test_that("rasterOverlap catches mismatched geometries", {
  r1 <- get_mock_raster(extent = c(0, 10, 0, 10))

  r2_bad_ext <- get_mock_raster(extent = c(5, 15, 5, 15))
  expect_error(rasterOverlap(r1, r2_bad_ext), "Rasters do not have the same geometry")

  r2_bad_res <- get_mock_raster(res = 0.5)
  expect_error(rasterOverlap(r1, r2_bad_res), "Rasters do not have the same geometry")
})

test_that("rasterOverlap handles negative values (log scale) with a warning", {
  r_pos <- get_mock_raster(vals = c(1, 2, 3, 4))
  r_neg <- get_mock_raster(vals = c(-1, -2, -3, -4))

  expect_warning(res1 <- rasterOverlap(r_neg, r_pos), "Negative values found in r1. Assuming log-scale")
  expect_warning(res2 <- rasterOverlap(r_pos, r_neg), "Negative values found in r2. Assuming log-scale")

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))
})

test_that("rasterOverlap catches rasters that cannot be normalized", {
  r1 <- get_mock_raster(vals = 1)
  r_zero <- get_mock_raster(vals = 0)

  expect_error(rasterOverlap(r1, r_zero), "One or more layers in r1 or r2 sum to 0 or NA")

  r_na <- get_mock_raster(vals = NA)

  expect_error(rasterOverlap(r1, r_na), "One or more layers in r1 or r2 sum to 0 or NA")
})

test_that("rasterOverlap correctly calculates affinity for identical distributions", {
  r1 <- get_mock_raster(vals = runif(100, 0.1, 1))

  affinity <- rasterOverlap(r1, r1)

  expect_equal(as.numeric(affinity), 1, tolerance = 1e-6)
})

test_that("rasterOverlap correctly calculates affinity for entirely disjoint distributions", {
  vals1 <- c(rep(1, 50), rep(0, 50))
  vals2 <- c(rep(0, 50), rep(1, 50))

  r1 <- get_mock_raster(vals = vals1)
  r2 <- get_mock_raster(vals = vals2)

  affinity <- rasterOverlap(r1, r2)

  expect_equal(as.numeric(affinity), 0, tolerance = 1e-6)
})

get_mock_getUD_output <- function(with_uncertainty = TRUE, log = FALSE, layers = 1) {
  # Base UD
  ud_list <- lapply(1:layers, function(i) {
    r <- get_mock_raster(vals = runif(100, 0.1, 1))
    if (log) terra::values(r) <- log(terra::values(r))
    return(r)
  })
  r_ud <- do.call(c, ud_list)
  names(r_ud) <- rep(ifelse(log, "log_UD", "UD"), layers)

  if (with_uncertainty) {
    # Add dummy SE and CV layers
    se_list <- lapply(1:layers, function(i) get_mock_raster(vals = runif(100, 0.01, 0.1)))
    r_se <- do.call(c, se_list)
    names(r_se) <- rep(ifelse(log, "log_UD_SE_delta", "UD_SE_delta"), layers)

    cv_list <- lapply(1:layers, function(i) get_mock_raster(vals = runif(100, 0.05, 0.2)))
    r_cv <- do.call(c, cv_list)
    names(r_cv) <- rep(ifelse(log, "log_UD_CV_delta", "UD_CV_delta"), layers)

    return(c(r_ud, r_se, r_cv))
  }

  return(r_ud)
}

test_that("rasterOverlap automatically filters out getUD uncertainty layers", {
  # r1 has 3 layers (UD, SE, CV), r2 has 1 layer (UD)
  r1 <- get_mock_getUD_output(with_uncertainty = TRUE, log = FALSE)
  r2 <- get_mock_getUD_output(with_uncertainty = FALSE, log = FALSE)

  # The smart filter should cleanly reduce r1 to just the UD layer
  expect_silent(res <- rasterOverlap(r1, r2))

  expect_length(res, 1)
  expect_null(names(res)) # Expecting unnamed vector for single-layer
  expect_true(res >= 0 && res <= 1)
})

test_that("rasterOverlap smart filter works correctly with log_UD", {
  r1 <- get_mock_getUD_output(with_uncertainty = TRUE, log = TRUE)
  r2 <- get_mock_getUD_output(with_uncertainty = FALSE, log = TRUE)

  # Catch BOTH warnings to prevent the second one from leaking into the console
  expect_warning(
    expect_warning(
      res <- rasterOverlap(r1, r2),
      "Negative values found in r1"
    ),
    "Negative values found in r2"
  )

  expect_length(res, 1)
  expect_null(names(res)) # Expecting unnamed vector for single-layer
  expect_true(res >= 0 && res <= 1)
})

test_that("rasterOverlap calculates pairwise affinities for multi-layer (dynamic) getUD outputs", {
  # Both are dynamic UDs with 2 time steps (6 layers total if with_uncertainty)
  r1 <- get_mock_getUD_output(with_uncertainty = TRUE, log = FALSE, layers = 2)
  r2 <- get_mock_getUD_output(with_uncertainty = FALSE, log = FALSE, layers = 2)

  # Smart filter should reduce r1 from 6 layers to 2 layers, matching r2
  res <- rasterOverlap(r1, r2)

  expect_length(res, 2)
  expect_named(res, c("UD", "UD")) # Expecting named vector for multi-layer
})

test_that("rasterOverlap calculates pairwise affinities for generic multi-layer rasters", {
  r1 <- c(get_mock_raster(vals = runif(100, 0.1, 1)), get_mock_raster(vals = runif(100, 0.1, 1)))
  names(r1) <- c("Habitat_A", "Habitat_B")

  r2 <- c(get_mock_raster(vals = runif(100, 0.1, 1)), get_mock_raster(vals = runif(100, 0.1, 1)))
  names(r2) <- c("Habitat_A", "Habitat_B")

  # Because neither has "UD" or "log_UD" in the name, the filter is bypassed
  res <- rasterOverlap(r1, r2)

  expect_length(res, 2)
  expect_named(res, c("Habitat_A", "Habitat_B"))
})

test_that("rasterOverlap catches mismatched layer counts AFTER smart filtering (User Trap)", {
  # r1 is a dynamic UD (2 time steps)
  r1 <- get_mock_getUD_output(with_uncertainty = FALSE, log = FALSE, layers = 2)

  # r2 is a static UD (1 time step)
  r2 <- get_mock_getUD_output(with_uncertainty = TRUE, log = FALSE, layers = 1)

  # After filtering, r1 has 2 layers and r2 has 1 layer. This should throw a clear error.
  expect_error(
    rasterOverlap(r1, r2),
    "Rasters must have the same number of layers"
  )
})

test_that("rasterOverlap catches mismatched layer counts for generic rasters (User Trap)", {
  r1 <- c(get_mock_raster(), get_mock_raster())
  names(r1) <- c("A", "B")

  r2 <- get_mock_raster()
  names(r2) <- "A"

  expect_error(
    rasterOverlap(r1, r2),
    "Rasters must have the same number of layers"
  )
})

test_that("rasterOverlap local argument correctly restricts affinity calculation", {
  r1 <- get_mock_raster(vals = 0.1)
  r2 <- get_mock_raster(vals = 0.9)

  # Make the rasters identical only in the center of the map
  center_ext <- terra::ext(4, 6, 4, 6)
  r1[center_ext] <- 1
  r2[center_ext] <- 1

  # Global affinity should be heavily penalized by the edges (< 1)
  global_affinity <- as.numeric(rasterOverlap(r1, r2))
  expect_true(global_affinity < 1)

  # Local affinity (cropped to the identical center) should be exactly 1
  local_affinity <- as.numeric(rasterOverlap(r1, r2, local = center_ext))
  expect_equal(local_affinity, 1, tolerance = 1e-6)
})

test_that("rasterOverlap local argument accepts SpatVector objects", {
  r1 <- get_mock_raster(vals = 0.1)
  r2 <- get_mock_raster(vals = 0.9)
  center_ext <- terra::ext(4, 6, 4, 6)
  r1[center_ext] <- 1
  r2[center_ext] <- 1

  # Convert the extent to a SpatVector polygon
  center_vec <- terra::as.polygons(center_ext)

  local_affinity <- as.numeric(rasterOverlap(r1, r2, local = center_vec))
  expect_equal(local_affinity, 1, tolerance = 1e-6)
})

test_that("rasterOverlap local argument accepts data frames and uses 'coord' correctly", {
  r1 <- get_mock_raster(vals = 0.1)
  r2 <- get_mock_raster(vals = 0.9)
  center_ext <- terra::ext(4, 6, 4, 6)
  r1[center_ext] <- 1
  r2[center_ext] <- 1

  # 1. Standard data frame with default coord = c("x", "y")
  df_standard <- data.frame(x = c(4, 6), y = c(4, 6))
  affinity_df <- as.numeric(rasterOverlap(r1, r2, local = df_standard))
  expect_equal(affinity_df, 1, tolerance = 1e-6)

  # 2. Simulated data frame with custom coord = c("mu.x", "mu.y")
  df_custom <- data.frame(mu.x = c(4, 6), mu.y = c(4, 6))
  affinity_custom <- as.numeric(rasterOverlap(r1, r2, local = df_custom, coord = c("mu.x", "mu.y")))
  expect_equal(affinity_custom, 1, tolerance = 1e-6)
})

test_that("rasterOverlap local argument automatically extracts estimated tracks from fitLangevin objects", {
  r1 <- get_mock_raster(vals = 0.1)
  r2 <- get_mock_raster(vals = 0.9)
  center_ext <- terra::ext(4, 6, 4, 6)
  r1[center_ext] <- 1
  r2[center_ext] <- 1

  # Construct a lightweight mock fitLangevin object
  mock_fit <- list(
    estimates = list(
      random = list(
        mu = list(
          est = data.frame(mu.x = c(4, 6), mu.y = c(4, 6))
        )
      )
    )
  )
  class(mock_fit) <- "fitLangevin"

  # It should bypass the coord argument and natively extract mu.x/mu.y
  affinity_fit <- as.numeric(rasterOverlap(r1, r2, local = mock_fit))
  expect_equal(affinity_fit, 1, tolerance = 1e-6)
})

test_that("rasterOverlap local argument throws informative errors for malformed inputs", {
  r1 <- get_mock_raster()
  r2 <- get_mock_raster()

  # 1. Data frame missing specified coord columns
  bad_df <- data.frame(lon = c(4, 6), lat = c(4, 6))
  expect_error(
    rasterOverlap(r1, r2, local = bad_df),
    "If 'local' is a data frame, it must contain the columns specified in 'coord' \\(currently: 'x' and 'y'\\)."
  )

  # 2. fitLangevin object lacking estimated random effects (mu)
  bad_fit <- list(estimates = list(natural = data.frame()))
  class(bad_fit) <- "fitLangevin"
  expect_error(
    rasterOverlap(r1, r2, local = bad_fit),
    "The provided 'fitLangevin' object does not contain estimated locations \\(mu\\)."
  )

  # 3. Completely unsupported object class
  expect_error(
    rasterOverlap(r1, r2, local = "just a string"),
    "'local' must be a SpatExtent, SpatVector, fitLangevin object, or a coordinate data frame."
  )
})
