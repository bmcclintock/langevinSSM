# tests/testthat/test_residuals.R

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

  # Fit the model so we have a baseline to compare against
  fit <- suppressMessages(suppressWarnings({
    fitLangevin(
      data = dat,
      model = "overdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      silent = TRUE
    )
  }))

  fit$residuals <- suppressMessages(residuals(fit, dat, spatialCovs))

  return(list(fit = fit, dat = dat, covs = spatialCovs))
}

# Generate the mock objects once for the test file
mock_env <- get_mock_res_fit()
fit_base <- mock_env$fit
dat_base <- mock_env$dat
covs_base <- mock_env$covs

test_that("residuals catches invalid inputs and legacy fit objects", {

  # 1. Non-fitLangevin object
  # To test our internal error handling, it MUST have the class to trigger S3 dispatch,
  # but lack the internal architecture (like tmb_setup).
  bad_fit <- list(par = c(1,2,3))
  class(bad_fit) <- "fitLangevin"

  expect_error(
    residuals(bad_fit, dat_base, covs_base),
    "does not contain a 'tmb_setup' blueprint"
  )

  # 2. Legacy fit object (missing the tmb_setup blueprint)
  legacy_fit <- fit_base
  legacy_fit$tmb_setup <- NULL
  expect_error(
    residuals(legacy_fit, dat_base, covs_base),
    "does not contain a 'tmb_setup' blueprint"
  )
})

test_that("residuals returns correctly formatted output", {

  # Run a fresh post-hoc extraction
  res_fresh <- suppressMessages(residuals(fit_base, dat_base, covs_base))

  # Check class
  expect_s3_class(res_fresh, "resLangevin")
  expect_s3_class(res_fresh, "data.frame")

  # Check required columns
  expect_true(all(c("id", "date", "residual.x", "residual.y") %in% names(res_fresh)))

  # Check rows perfectly match the input data
  expect_equal(nrow(res_fresh), nrow(dat_base))
})

test_that("residuals successfully accepts and passes alternative OSA methods", {

  # Reconstruct using "fullGaussian" instead of the default
  res_full <- suppressMessages(
    residuals(fit_base, dat_base, covs_base, method = "fullGaussian")
  )

  expect_s3_class(res_full, "resLangevin")

  # The method should run without error and produce valid numeric residuals.
  # Note: For this specific Gaussian model, fullGaussian and oneStepGaussianOffMode
  # are mathematically identical, so we just verify it produced valid numerics.
  valid_residuals <- na.omit(res_full$residual.x)
  expect_true(length(valid_residuals) > 0)
  expect_true(is.numeric(valid_residuals))
})

test_that("residuals handles individual track failures gracefully without crashing the whole process", {

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
    res_sabotaged <- suppressMessages(residuals(bad_fit, bad_dat, covs_base)),
    "OSA calculation failed for track ID 2\\. Error: Not enough valid observations"
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

test_that("underdamped fitLangevin calculates OSA residuals correctly with duplicated times, NA gaps, and known locations", {

  # 1. Simulate Dummy Spatial Covariate
  set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-500, 500, -500, 500))
  terra::values(r) <- runif(terra::ncell(r))
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2. Simulate Track Data using simLangevin
  init_par <- list(beta = 0, sigma = 5, gamma=0.5)

  dat <- suppressMessages(simLangevin(obsPerAnimal=100,
                                      model = "underdamped",
                                      spatialCovs = spatialCovs,
                                      par = init_par
  ))

  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dat$date <- start_time + (1:nrow(dat)) * 3600

  # FORCEFULLY assign observation errors to everything initially
  dat$x.err <- 2
  dat$y.err <- 2
  dat$lc <- "G"
  dat$id <- "A"

  # 3. Inject Edge Cases

  # A. Hard Duplicate 1 (With observation error)
  dup_row1 <- dat[5, ]
  dup_row1$x <- dup_row1$x + 1.0 # Slight GPS jitter
  dat <- rbind(dat, dup_row1)

  # B. Hard Duplicate 2 (Known location, NO observation error)
  # We set x.err and y.err to NA to trigger the "perfect track" logic.
  dat$x.err[20] <- NA
  dat$y.err[20] <- NA
  dup_row2 <- dat[20, ]
  dat <- rbind(dat, dup_row2)

  dat <- dat[order(dat$date), ]

  # C. True Interpolation Point (NA coordinates)
  pt <- data.frame(id = dat$id[1], date = dat$date[10] + 1800)

  class(dat) <- "data.frame"
  fmt_dat <- suppressWarnings(formatData(dat, predTimes = pt, time.unit = "hours"))

  # 4. Fit Model
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "underdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      silent = TRUE
    )
  }))

  fit$residuals <- suppressMessages(residuals(fit, fmt_dat, spatialCovs))

  # 5. Assertions
  expect_s3_class(fit, "fitLangevin")
  expect_true(is.data.frame(fit$residuals))

  # Total rows should match formatted data exactly
  expect_equal(nrow(fit$residuals), nrow(fmt_dat))

  # 1. The interpolation point (NA coordinates) should have NA residuals
  na_obs_idx <- which(is.na(fmt_dat$x))
  expect_true(length(na_obs_idx) == 1)
  expect_true(all(is.na(fit$residuals$residual.x[na_obs_idx])))

  # 2. The first observation of EVERY track should have NA residuals
  first_obs_idx <- which(!duplicated(fmt_dat$id))
  expect_true(all(is.na(fit$residuals$residual.x[first_obs_idx])))

  # 3. Known Locations (No observation error)
  # Because these have no error, TMB skips their observation likelihood (isd=0).
  # With no likelihood to evaluate, their residuals MUST be NA.
  known_obs_idx <- which(is.na(fmt_dat$x.err) & !is.na(fmt_dat$x))
  expect_true(length(known_obs_idx) == 2) # The original + the duplicate
  expect_true(all(is.na(fit$residuals$residual.x[known_obs_idx])))

  # 4. ALL OTHER real observations
  # This includes the hard duplicate WITH observation error! Since it has a valid
  # observation variance, TMB successfully computes a numeric residual for it.
  real_obs_idx <- which(!is.na(fmt_dat$x))

  # Filter out the first observation and the known (error-free) observations
  obs_to_check <- setdiff(real_obs_idx, c(first_obs_idx, known_obs_idx))

  expect_true(all(!is.na(fit$residuals$residual.x[obs_to_check])))
  expect_true(all(!is.na(fit$residuals$residual.y[obs_to_check])))
})

test_that("overdamped fitLangevin calculates OSA residuals correctly with duplicated times, NA gaps, and known locations", {

  # 1. Simulate Dummy Spatial Covariate
  set.seed(10,kind="Mersenne-Twister",normal.kind="Inversion")
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-500, 500, -500, 500))
  terra::values(r) <- runif(terra::ncell(r))
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2. Simulate Track Data using simLangevin
  init_par <- list(beta = 0, sigma = 5)

  dat <- suppressMessages(simLangevin(obsPerAnimal=100,
                                      model = "overdamped",
                                      spatialCovs = spatialCovs,
                                      par = init_par
  ))

  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dat$date <- start_time + (1:nrow(dat)) * 3600

  # FORCEFULLY assign observation errors to everything initially
  dat$x.err <- 2
  dat$y.err <- 2
  dat$lc <- "G"
  dat$id <- "A"

  # 3. Inject Edge Cases

  # A. Hard Duplicate 1 (With observation error)
  dup_row1 <- dat[5, ]
  dup_row1$x <- dup_row1$x + 1.0 # Slight GPS jitter
  dat <- rbind(dat, dup_row1)

  # B. Hard Duplicate 2 (Known location, NO observation error)
  # We set x.err and y.err to NA to trigger the "perfect track" logic.
  dat$x.err[20] <- NA
  dat$y.err[20] <- NA
  dup_row2 <- dat[20, ]
  dat <- rbind(dat, dup_row2)

  dat <- dat[order(dat$date), ]

  # C. True Interpolation Point (NA coordinates)
  pt <- data.frame(id = dat$id[1], date = dat$date[10] + 1800)

  class(dat) <- "data.frame"
  fmt_dat <- suppressWarnings(formatData(dat, predTimes = pt, time.unit = "hours"))

  # 4. Fit Model
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "overdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      silent = TRUE
    )
  }))

  fit$residuals <- suppressMessages(residuals(fit, fmt_dat, spatialCovs))

  # 5. Assertions
  expect_s3_class(fit, "fitLangevin")
  expect_true(is.data.frame(fit$residuals))

  # Total rows should match formatted data exactly
  expect_equal(nrow(fit$residuals), nrow(fmt_dat))

  # 1. The interpolation point (NA coordinates) should have NA residuals
  na_obs_idx <- which(is.na(fmt_dat$x))
  expect_true(length(na_obs_idx) == 1)
  expect_true(all(is.na(fit$residuals$residual.x[na_obs_idx])))

  # 2. Known Locations (No observation error) MUST be NA
  known_obs_idx <- which(is.na(fmt_dat$x.err) & !is.na(fmt_dat$x))
  expect_true(length(known_obs_idx) == 2) # The original + the duplicate
  expect_true(all(is.na(fit$residuals$residual.x[known_obs_idx])))

  # 3. The first observation of EVERY track MUST be NA
  # (Because we explicitly conditioned on it in fitLangevin.R!)
  first_obs_idx <- which(!duplicated(fmt_dat$id))
  expect_true(all(is.na(fit$residuals$residual.x[first_obs_idx])))
  expect_true(all(is.na(fit$residuals$residual.y[first_obs_idx])))

  # 4. ALL OTHER real observations (including the hard duplicate WITH error)
  real_obs_idx <- which(!is.na(fmt_dat$x))

  # Filter out the known (error-free) observations and the first observation
  obs_to_check <- setdiff(real_obs_idx, c(first_obs_idx, known_obs_idx))

  expect_true(all(!is.na(fit$residuals$residual.x[obs_to_check])))
  expect_true(all(!is.na(fit$residuals$residual.y[obs_to_check])))
})

test_that("print.resLangevin outputs summary and GoF tests correctly", {

  # Run a fresh extraction to ensure tests attribute is attached
  res_obj <- suppressMessages(residuals(fit_base, dat_base, covs_base))

  # Tiny mock datasets sometimes fail GoF tests and return NULL.
  # Inject a dummy test dataframe if necessary just to ensure the print method logic fires.
  if(is.null(attr(res_obj, "tests"))) {
    attr(res_obj, "tests") <- data.frame(ks_test_x = 0.99)
  }

  # Capture the console output of the print method explicitly to bypass testthat scoping
  printed_output <- capture.output(print.resLangevin(res_obj))

  # Verify Header
  expect_true(any(grepl("=== One-Step-Ahead \\(OSA\\) Residuals ===", printed_output)))

  # Verify Structural Info
  expect_true(any(grepl("Total observations:", printed_output)))
  expect_true(any(grepl("Number of tracks:", printed_output)))

  # Verify GoF Test Table was printed
  expect_true(any(grepl("--- Goodness-of-Fit Tests ---", printed_output)))
  expect_true(any(grepl("KS_x", printed_output)))

  # Verify Residual Summary was printed
  expect_true(any(grepl("Residual Summary:", printed_output)))
  expect_true(any(grepl("residual.x", printed_output)))

  # Verify the plotting hint is present
  #expect_true(any(grepl("Tip: Use plot\\(\\)", printed_output)))
})
