# tests/testthat/test_plotRaster.R

test_that("plotRaster returns a ggplot object", {
  r1 <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10), vals = runif(100))
  names(r1) <- "habitat"

  p <- plotRaster(r1)
  expect_s3_class(p, "ggplot")
})

test_that("plotRaster correctly applies custom legend titles", {
  r1 <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10), vals = runif(100))

  p <- plotRaster(r1, legend.title = "Custom Scale")
  expect_s3_class(p, "ggplot")

  # Dig into ggplot internals to verify the scale name was set
  scale_names <- sapply(p$scales$scales, function(x) x$name)
  expect_true("Custom Scale" %in% scale_names)
})
