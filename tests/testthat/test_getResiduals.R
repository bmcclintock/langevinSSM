# tests/testthat/test_getResiduals.R

# --- Helper to generate a fast, reusable mock fit ---
get_mock_res_fit <- function() {
  set.seed(42)
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- runif(100)
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2 animals, 15 obs each to test multi-track looping
  init_par <- list(beta = 0, sigma = 5)
  dat <- suppressMessages(simLangevin(
    obsPerAnimal = 15,
    nbAnimals = 2,
    model = "overdamped",
    spatialCovs = spatialCovs,
    par = init_par,
    measurementError = list(x.sd = 2, y.sd = 2)
  ))

  # Ensure valid dates
  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dat$date <- start_time + (1:nrow(dat)) * 3600

  # Fit the model with calcResiduals = TRUE so we have a baseline to compare against
  fit <- suppressMessages(suppressWarnings({
    fitLangevin(
      data = dat,
      model = "overdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      calcResiduals = TRUE,
      silent = TRUE
    )
  }))

  return(list(fit = fit, dat = dat, covs = spatialCovs))
}

# Generate the mock objects once for the test file
mock_env <- get_mock_res_fit()
fit_base <- mock_env$fit
dat_base <- mock_env$dat
covs_base <- mock_env$covs

test_that("getResiduals catches invalid inputs and legacy fit objects", {

  # 1. Non-fitLangevin object
  bad_fit <- list(par = c(1,2,3))
  expect_error(
    getResiduals(bad_fit, dat_base, covs_base),
    "'fit' must be a fitLangevin object."
  )

  # 2. Legacy fit object (missing the tmb_setup blueprint)
  legacy_fit <- fit_base
  legacy_fit$tmb_setup <- NULL
  expect_error(
    getResiduals(legacy_fit, dat_base, covs_base),
    "does not contain a 'tmb_setup' blueprint"
  )
})

test_that("getResiduals returns correctly formatted output", {

  # Run a fresh post-hoc extraction
  res_fresh <- suppressMessages(getResiduals(fit_base, dat_base, covs_base))

  # Check class
  expect_s3_class(res_fresh, "resLangevin")
  expect_s3_class(res_fresh, "data.frame")

  # Check required columns
  expect_true(all(c("id", "date", "residual.x", "residual.y") %in% names(res_fresh)))

  # Check rows perfectly match the input data
  expect_equal(nrow(res_fresh), nrow(dat_base))
})

test_that("getResiduals post-hoc reconstruction perfectly matches fitLangevin internal calc", {

  # fit_base$residuals was calculated INSIDE fitLangevin using the active TMB obj.
  # res_fresh is calculated by RECONSTRUCTING the TMB obj from the blueprint.
  res_fresh <- suppressMessages(getResiduals(fit_base, dat_base, covs_base))

  # They should be mathematically identical
  expect_equal(res_fresh$residual.x, fit_base$residuals$residual.x)
  expect_equal(res_fresh$residual.y, fit_base$residuals$residual.y)
})

test_that("getResiduals successfully accepts and passes alternative OSA methods", {

  # Reconstruct using "fullGaussian" instead of the default
  res_full <- suppressMessages(
    getResiduals(fit_base, dat_base, covs_base, method = "fullGaussian")
  )

  expect_s3_class(res_full, "resLangevin")

  # The method should run without error and produce valid numeric residuals.
  # Note: For this specific Gaussian model, fullGaussian and oneStepGaussianOffMode
  # are mathematically identical, so we just verify it produced valid numerics.
  valid_residuals <- na.omit(res_full$residual.x)
  expect_true(length(valid_residuals) > 0)
  expect_true(is.numeric(valid_residuals))
})

test_that("getResiduals handles individual track failures gracefully without crashing the whole process", {

  bad_dat <- dat_base

  # Sabotage Track 2 by leaving only its FIRST observation intact.
  # This triggers the "not enough valid observations" warning, but anchors
  # the random walk so TMB's precision matrix doesn't become singular and crash Track 1.
  t2_idx <- which(bad_dat$id == "2")
  bad_dat$x[t2_idx[-1]] <- NA
  bad_dat$y[t2_idx[-1]] <- NA

  # need to update the fit object's signature so it allows this
  bad_fit <- fit_base
  bad_fit$signatures$data <- langevinSSM:::get_data_signature(bad_dat)

  # Expect a warning that Track 2 didn't have enough observations
  expect_warning(
    res_sabotaged <- suppressMessages(getResiduals(bad_fit, bad_dat, covs_base)),
    "Track ID 2 does not have enough valid observations"
  )

  # BUT the function should not have crashed. It should still return an object...
  expect_s3_class(res_sabotaged, "resLangevin")

  # ...where Animal 1 successfully got its residuals...
  anim_1_resids <- na.omit(res_sabotaged$residual.x[res_sabotaged$id == "1"])
  expect_true(length(anim_1_resids) > 0)

  # ...and Animal 2 is entirely NA
  anim_2_resids <- na.omit(res_sabotaged$residual.x[res_sabotaged$id == "2"])
  expect_equal(length(anim_2_resids), 0)
})
