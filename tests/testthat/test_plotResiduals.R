# tests/testthat/test_plotResiduals.R

# --- Helper: Generate a fast mock resLangevin object ---
get_mock_resLangevin <- function() {
  set.seed(42, kind="Mersenne-Twister", normal.kind="Inversion")

  df <- data.frame(
    id = rep(c("A", "B"), each = 50),
    date = as.POSIXct("2024-01-01 00:00:00", tz = "UTC") + (1:100) * 3600,
    residual.x = rnorm(100),
    residual.y = rnorm(100)
  )

  # Inject a few NAs (simulating the first obs of a track) to test robust na.omit handling
  df$residual.x[c(1, 51)] <- NA
  df$residual.y[c(1, 51)] <- NA

  class(df) <- c("resLangevin", "data.frame")
  return(df)
}

test_that("plotResiduals catches invalid inputs", {
  # Standard data frame without the resLangevin class
  bad_res <- data.frame(id = "A", residual.x = 1, residual.y = 1)

  expect_error(
    plotResiduals(bad_res),
    "must be a 'resLangevin' object"
  )
})

test_that("plotResiduals generates aggregated plots by default", {
  res <- get_mock_resLangevin()
  p_list <- plotResiduals(res)

  # Should return a single, flat list of 6 plots
  expect_type(p_list, "list")
  expect_equal(length(p_list), 6)
  expect_named(p_list, c("qq_x", "qq_y", "acf_x", "acf_y", "qq_mah", "acf_mah"))

  # Ensure they are actually ggplot objects
  expect_s3_class(p_list$qq_x, "ggplot")
  expect_s3_class(p_list$acf_mah, "ggplot")
})

test_that("plotResiduals generates separate nested plots when tracks='all'", {
  res <- get_mock_resLangevin()
  p_list <- plotResiduals(res, tracks = "all")

  # Should return a list of lists, named by ID
  expect_type(p_list, "list")
  expect_named(p_list, c("A", "B"))

  # Check inner structure
  expect_equal(length(p_list$A), 6)
  expect_s3_class(p_list$A$qq_x, "ggplot")
  expect_s3_class(p_list$B$acf_mah, "ggplot")
})

test_that("plotResiduals filters specific tracks correctly", {
  res <- get_mock_resLangevin()
  p_list <- plotResiduals(res, tracks = c("B"))

  # Should return a list containing ONLY track B
  expect_type(p_list, "list")
  expect_named(p_list, c("B"))
  expect_null(p_list$A)
})

test_that("plotResiduals handles missing/invalid track IDs gracefully", {
  res <- get_mock_resLangevin()

  # Requesting a track that doesn't exist should throw a warning, not crash
  expect_warning(
    p_list <- plotResiduals(res, tracks = c("C")),
    "No data found for track ID: C"
  )

  # Because "C" was the only request, the resulting list should be empty
  expect_equal(length(p_list), 0)

  # Requesting a mix of valid and invalid tracks
  expect_warning(
    p_list_mixed <- plotResiduals(res, tracks = c("A", "C")),
    "No data found for track ID: C"
  )

  # It should successfully plot A, and skip C
  expect_named(p_list_mixed, c("A"))
})
