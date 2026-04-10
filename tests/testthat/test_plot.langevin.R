# tests/testthat/test_plot.langevin.R

# --- Helper functions for mock data ---
get_mock_covs <- function(dynamic = FALSE) {
  r1 <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10), vals = runif(100))
  if (dynamic) {
    r1 <- c(r1, r1 * 0.5)
    terra::time(r1) <- as.POSIXct(c("2023-01-01 00:00:00", "2023-01-02 00:00:00"), tz = "UTC")
  }
  list(cov1 = r1)
}

# Generate a sequence of dates bridging the dynamic raster times
mock_dates <- rep(seq(
  from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
  to = as.POSIXct("2023-01-02 00:00:00", tz = "UTC"),
  length.out = 5
), 2)

get_mock_fit <- function() {
  fit <- list(
    estimates = list(
      natural = matrix(c(1.5), nrow = 1, dimnames = list(c("beta"), "Estimate")),
      random = list(
        mu = list(
          est = data.frame(
            mu.x = c(2.1, 2.6, 3.1, 3.6, 4.1, 6.1, 6.6, 7.1, 7.6, 8.1),
            mu.y = c(2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, 7.5, 8.0),
            id = rep(c("A", "B"), each = 5),
            date = mock_dates
          )
        )
      )
    ),
    signatures = list(
      data = NULL,
      covs = NULL
    )
  )
  class(fit) <- "fitLangevin"
  return(fit)
}

get_mock_data <- function() {
  # 5 points per track, now strictly including dates
  dat <- data.frame(
    id = rep(c("A", "B"), each = 5),
    date = mock_dates,
    dt = rep(6, 10),
    x = c(2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, 7.5, 8.0),
    y = c(2.1, 2.4, 3.1, 3.6, 4.1, 6.1, 6.4, 7.1, 7.6, 8.1),
    smaj = NA, smin = NA, eor = NA, x.sd = NA, y.sd = NA
  )
  # Actually route it through the real validator so we know the mock data is legal!
  dat <- class_dataLangevin(dat, time.unit = "hours")
  return(dat)
}

get_mock_sim <- function() {
  # 5 points per track, with dates and latent states
  dat <- data.frame(
    id = rep(c("A", "B"), each = 5),
    date = mock_dates,
    dt = rep(6, 10),
    x = c(2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, 7.5, 8.0),
    y = c(2.1, 2.4, 3.1, 3.6, 4.1, 6.1, 6.4, 7.1, 7.6, 8.1),
    mu.x = c(2.1, 2.6, 3.1, 3.6, 4.1, 6.1, 6.6, 7.1, 7.6, 8.1),
    mu.y = c(2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, 7.5, 8.0),
    smaj = NA, smin = NA, eor = NA, x.sd = NA, y.sd = NA
  )
  dat <- class_dataLangevin(dat, time.unit = "hours")
  class(dat) <- unique(c("simLangevin", class(dat)))
  return(dat)
}

# ---------------------------------------------------------
# Tests for plot.fitLangevin
# ---------------------------------------------------------

test_that("plot.fitLangevin catches missing spatialCovs", {
  fit <- get_mock_fit()
  expect_error(plot(fit), "You must provide the 'spatialCovs'")
})

test_that("plot.fitLangevin works for static rasters (Single Layer)", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)

  p <- plot(fit, spatialCovs = covs)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"), info = "Static raster should not use facet_wrap")
})

test_that("plot.fitLangevin works with observed data overlay", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)
  dat <- get_mock_data()

  p <- plot(fit, spatialCovs = covs, data = dat)
  expect_s3_class(p, "ggplot")
})

test_that("plot.fitLangevin handles dynamic rasters with and without time subsetting", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = TRUE)

  # All layers (faceting should trigger)
  p_all <- plot(fit, spatialCovs = covs)
  expect_s3_class(p_all, "ggplot")
  expect_true(inherits(p_all$facet, "FacetWrap"), info = "Multi-layer UD should trigger facet_wrap")

  # Subsetting by time index (reduces to single layer)
  p_sub <- plot(fit, spatialCovs = covs, time = 1)
  expect_s3_class(p_sub, "ggplot")
  expect_true(inherits(p_sub$facet, "FacetNull"), info = "Subsetting to a single time layer should remove facet_wrap")
})

test_that("plot.fitLangevin compact = FALSE returns a list of individual plots", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)

  p_list <- plot(fit, spatialCovs = covs, compact = FALSE)

  # Because the mock data has 2 IDs ("A" and "B"), it should return a list of 2 ggplots
  expect_true(is.list(p_list))
  expect_equal(length(p_list), 2)
  expect_s3_class(p_list[[1]], "ggplot")
  expect_s3_class(p_list[[2]], "ggplot")
})

# ---------------------------------------------------------
# Tests for plot.dataLangevin
# ---------------------------------------------------------

test_that("plot.dataLangevin catches missing spatialCovs", {
  dat <- get_mock_data()
  expect_error(plot(dat), "You must provide the 'spatialCovs'")
})

test_that("plot.dataLangevin returns a list of plots per covariate (Single Layer)", {
  dat <- get_mock_data()
  covs <- get_mock_covs(dynamic = FALSE)
  covs$cov2 <- covs$cov1 # Add a second covariate to test list length

  p_list <- plot(dat, spatialCovs = covs)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 2)
  expect_s3_class(p_list[["cov1"]], "ggplot")
  expect_s3_class(p_list[["cov2"]], "ggplot")
  expect_true(inherits(p_list[["cov1"]]$facet, "FacetNull"))
})

test_that("plot.dataLangevin compact = FALSE returns nested list of plots", {
  dat <- get_mock_data()
  covs <- get_mock_covs(dynamic = FALSE)

  p_nested <- plot(dat, spatialCovs = covs, compact = FALSE)

  expect_true(is.list(p_nested[["cov1"]]))
  expect_equal(length(p_nested[["cov1"]]), 2) # Two track IDs ("A" and "B")
  expect_s3_class(p_nested[["cov1"]][["A"]], "ggplot")
})

test_that("plot.dataLangevin handles multi-layer (dynamic) rasters", {
  dat <- get_mock_data()
  covs <- get_mock_covs(dynamic = TRUE)

  # Only 1 covariate provided, so it returns 1 ggplot object in the list.
  # But when plotted, that single ggplot will render with facet_wrap for the multiple layers.
  p_list <- plot(dat, spatialCovs = covs)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 1)
  expect_s3_class(p_list[["cov1"]], "ggplot")
  expect_true(inherits(p_list[["cov1"]]$facet, "FacetWrap"), info = "Multi-layer covariate should trigger facet_wrap")
})

test_that("plot.dataLangevin handles time subsetting with mixed static/dynamic covs without warning", {
  dat <- get_mock_data()
  covs <- get_mock_covs(dynamic = TRUE)
  covs$static_cov <- get_mock_covs(dynamic = FALSE)$cov1 # Mix of dynamic and static

  # Because at least one covariate is dynamic, NO warning should be produced!
  p_list <- plot(dat, spatialCovs = covs, time = 1)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 2)
  expect_true(inherits(p_list[["cov1"]]$facet, "FacetNull"), info = "Dynamic raster should be reduced to single layer")
  expect_true(inherits(p_list[["static_cov"]]$facet, "FacetNull"), info = "Static raster remains single layer without crashing")
})

test_that("plot.dataLangevin warns if time provided but ALL covariates are static", {
  dat <- get_mock_data()
  covs <- get_mock_covs(dynamic = FALSE) # Only static
  covs$static_cov2 <- get_mock_covs(dynamic = FALSE)$cov1

  expect_warning(
    plot(dat, spatialCovs = covs, time = 1),
    "all spatial covariates are static"
  )
})

# ---------------------------------------------------------
# Tests for plot.simLangevin
# ---------------------------------------------------------

test_that("plot.simLangevin without beta falls back to dataLangevin multi-covariate behavior", {
  sim <- get_mock_sim()
  covs <- get_mock_covs(dynamic = FALSE)
  covs$cov2 <- covs$cov1

  p_list <- plot(sim, spatialCovs = covs)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 2)
  expect_s3_class(p_list[["cov1"]], "ggplot")
})

test_that("plot.simLangevin with beta calculates UD and returns a single combined plot", {
  sim <- get_mock_sim()
  covs <- get_mock_covs(dynamic = FALSE)

  p <- plot(sim, spatialCovs = covs, beta = c(1.5))

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"))
})

test_that("plot.simLangevin handles multi-layer (dynamic) rasters without beta", {
  sim <- get_mock_sim()
  covs <- get_mock_covs(dynamic = TRUE)

  # Falling back to dataLangevin behavior
  p_list <- plot(sim, spatialCovs = covs)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 1)
  expect_s3_class(p_list[["cov1"]], "ggplot")
  expect_true(inherits(p_list[["cov1"]]$facet, "FacetWrap"))
})

test_that("plot.simLangevin handles multi-layer (dynamic) rasters with beta", {
  sim <- get_mock_sim()
  covs <- get_mock_covs(dynamic = TRUE)

  # Calculates a dynamic UD, returning one faceted ggplot
  p <- plot(sim, spatialCovs = covs, beta = c(1.5))

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"), info = "Multi-layer UD calculation should trigger facet_wrap")
})

test_that("plot.simLangevin handles time subsetting with beta", {
  sim <- get_mock_sim()
  covs <- get_mock_covs(dynamic = TRUE)

  p <- plot(sim, spatialCovs = covs, beta = c(1.5), time = 1)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"), info = "Dynamic UD should be reduced to single layer")
})

# ---------------------------------------------------------
# Tests for plot.udLangevin
# ---------------------------------------------------------

# --- Helper function for mock udLangevin data ---
get_mock_ud <- function(dynamic = FALSE, with_uncertainty = FALSE) {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10), vals = runif(100))
  names(r) <- "UD"

  if (dynamic) {
    r <- c(r, r * 0.5)
    terra::time(r) <- as.POSIXct(c("2023-01-01 00:00:00", "2023-01-02 00:00:00"), tz = "UTC")
  }

  out <- list(UD = r)

  if (with_uncertainty) {
    # Create arbitrary SE and CV rasters based on the UD
    se <- r * 0.1
    names(se) <- rep("UD_SE", terra::nlyr(se))
    cv <- r * 0.05
    names(cv) <- rep("UD_CV", terra::nlyr(cv))

    out$SE <- se
    out$CV <- cv
  }

  # Using the class constructor format
  class(out) <- unique(c("udLangevin", class(out)))
  return(out)
}

test_that("plot.udLangevin works for basic UD without uncertainty", {
  ud <- get_mock_ud(with_uncertainty = FALSE)
  p <- plot(ud)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"), info = "Static base UD should not use facet_wrap")
  expect_equal(p$labels$title, "Utilization Distribution (UD)")
})

test_that("plot.udLangevin returns list of plots when SE and CV are present", {
  ud <- get_mock_ud(with_uncertainty = TRUE)
  p_list <- plot(ud)

  expect_type(p_list, "list")
  expect_equal(length(p_list), 3)
  expect_named(p_list, c("UD", "SE", "CV"))

  expect_s3_class(p_list[["UD"]], "ggplot")
  expect_s3_class(p_list[["SE"]], "ggplot")
  expect_s3_class(p_list[["CV"]], "ggplot")

  expect_equal(p_list[["SE"]]$labels$title, "UD log standard error (SE)")
  expect_equal(p_list[["CV"]]$labels$title, "UD coefficient of variation (CV)")
})

test_that("plot.udLangevin handles log = FALSE for SE", {
  ud <- get_mock_ud(with_uncertainty = TRUE)
  p_list <- plot(ud, log = FALSE)

  expect_type(p_list, "list")
  expect_s3_class(p_list[["SE"]], "ggplot")

  # Check that the title changed to indicate logarithmic scaling
  expect_equal(p_list[["SE"]]$labels$title, "UD standard error (SE)")
})

test_that("plot.udLangevin handles dynamic rasters with and without time subsetting", {
  ud <- get_mock_ud(dynamic = TRUE, with_uncertainty = TRUE)

  # All layers (faceting should trigger on all generated plots)
  p_all <- plot(ud)
  expect_true(inherits(p_all[["UD"]]$facet, "FacetWrap"), info = "Multi-layer UD should trigger facet_wrap")
  expect_true(inherits(p_all[["SE"]]$facet, "FacetWrap"), info = "Multi-layer SE should trigger facet_wrap")
  expect_true(inherits(p_all[["CV"]]$facet, "FacetWrap"), info = "Multi-layer CV should trigger facet_wrap")

  # Subsetting by time index (reduces to single layer)
  p_sub <- plot(ud, time = 1)
  expect_true(inherits(p_sub[["UD"]]$facet, "FacetNull"), info = "Subsetting dynamic UD to a single time layer should remove facet_wrap")
  expect_true(inherits(p_sub[["SE"]]$facet, "FacetNull"), info = "Subsetting dynamic SE to a single time layer should remove facet_wrap")
})
