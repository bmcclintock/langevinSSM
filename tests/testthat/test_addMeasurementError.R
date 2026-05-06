# tests/testthat/test_addMeasurementError.R

# --- Helpers ---
get_base_data <- function() {
  df <- data.frame(
    id = rep("A", 10),
    date = as.POSIXct("2024-01-01 00:00:00", tz = "UTC") + (1:10) * 3600,
    mu.x = seq(0, 100, length.out = 10),
    mu.y = seq(0, 100, length.out = 10)
  )
  class(df) <- c("dataLangevin", class(df))
  return(df)
}

# --- Tests: User Error Catching ---

test_that("addMeasurementError catches missing inputs", {
  df <- get_base_data()

  # Neither known error columns nor measurementError list provided
  expect_error(addMeasurementError(df), "No measurement error parameters provided")
})

test_that("addMeasurementError catches parameter scale conflicts", {
  df <- get_base_data()
  err <- list(x.sd = 2, y.sd = 2)

  # Conflicting psi
  expect_error(addMeasurementError(df, par = list(psi = 1, l_psi = 0), measurementError = err),
               "Cannot provide both 'psi' and 'l_psi'")

  # Conflicting tau
  expect_error(addMeasurementError(df, par = list(tau = c(1,1), l_tau = c(0,0)), measurementError = err),
               "Cannot provide both 'tau' and 'l_tau'")

  # Conflicting rho_o
  expect_error(addMeasurementError(df, par = list(rho_o = 0, l_rho_o = 0), measurementError = err),
               "Cannot provide both 'rho_o' and 'l_rho_o'")
})

test_that("addMeasurementError catches conflicting error types in measurementError list", {
  df <- get_base_data()

  # Providing both KF and LS generating parameters
  bad_err <- list(smaj.sd = 1, smin.sd = 0.5, x.sd = 1, y.sd = 1)

  expect_error(addMeasurementError(df, measurementError = bad_err),
               "Cannot provide both error ellipse")

  # Not a list
  expect_error(addMeasurementError(df, measurementError = c(x.sd = 1, y.sd = 1)),
               "must be a list")
})


# --- Tests: Generating New Errors ---

test_that("addMeasurementError generates new KF errors correctly", {
  df <- get_base_data()
  err <- list(smaj.sd = 2, smin.sd = 1, eor.lim = c(0, 180))

  set.seed(123)
  res <- suppressWarnings(addMeasurementError(df, measurementError = err))

  # Check columns created
  expect_true(all(c("x", "y", "smaj", "smin", "eor", "x.err", "y.err") %in% names(res)))

  # Check data types and structure
  expect_true(is.numeric(res$smaj))
  expect_true(all(is.na(res$x.err)))
  expect_true(all(res$lc == "G"))

  # Check noise was actually added (x should not perfectly equal mu.x anymore)
  expect_false(isTRUE(all.equal(res$x, res$mu.x)))

  # eor should be bounded between 0 and pi (since input eor.lim was in degrees, rcpp handles the rad conversion)
  expect_true(all(res$eor >= 0 & res$eor <= pi, na.rm = TRUE))
})

test_that("addMeasurementError generates new LS errors correctly", {
  df <- get_base_data()
  err <- list(x.sd = 5, y.sd = 5)

  set.seed(123)
  res <- suppressWarnings(addMeasurementError(df, measurementError = err))

  # Check columns
  expect_true(all(!is.na(res$x.err) & !is.na(res$y.err)))
  expect_true(all(is.na(res$smaj) & is.na(res$smin) & is.na(res$eor)))

  # Check noise was added
  expect_false(isTRUE(all.equal(res$y, res$mu.y)))
})


# --- Tests: Applying Known Errors ---

test_that("addMeasurementError applies known KF errors correctly", {
  df <- get_base_data()

  # Pre-populate errors directly into the dataframe
  df$smaj <- 5
  df$smin <- 2
  df$eor <- 1.5

  set.seed(456)
  res <- suppressWarnings(addMeasurementError(df)) # measurementError list intentionally omitted

  # x and y should be generated based on the hardcoded smaj/smin
  expect_false(isTRUE(all.equal(res$x, res$mu.x)))
  expect_true(all(!is.na(res$x)))

  # Original error columns should remain completely untouched
  expect_true(all(res$smaj == 5))
  expect_true(all(res$smin == 2))

  # LS should be NA
  expect_true(all(is.na(res$x.err)))
})

test_that("addMeasurementError applies known LS errors correctly", {
  df <- get_base_data()

  # Pre-populate errors directly
  df$x.err <- 5
  df$y.err <- 5

  set.seed(456)
  res <- suppressWarnings(addMeasurementError(df))

  # Noise added
  expect_false(isTRUE(all.equal(res$y, res$mu.y)))
  expect_true(all(!is.na(res$y)))

  # Original sd columns remain intact
  expect_true(all(res$x.err == 5))

  # KF columns should be NA
  expect_true(all(is.na(res$smin)))
})


# --- Tests: Parameter Parsing Verification ---

test_that("addMeasurementError properly parses natural and working parameter scales", {
  df <- get_base_data()
  df$x.err <- 1
  df$y.err <- 1

  # Using natural scale (tau = 10 scales the standard deviation x10)
  set.seed(111)
  res_natural <- suppressWarnings(addMeasurementError(df, par = list(tau = c(10, 10))))

  # Using working scale (l_tau = log(10))
  set.seed(111)
  res_working <- suppressWarnings(addMeasurementError(df, par = list(l_tau = c(log(10), log(10)))))

  # The mathematical outcome of the Rcpp multivariate normal generation should be identical
  expect_equal(res_natural$x, res_working$x)
  expect_equal(res_natural$y, res_working$y)
})

test_that("addMeasurementError catches conflicting error types in data (known errors)", {
  df <- get_base_data()

  # Pre-populate BOTH KF and LS error columns
  df$smaj <- 5
  df$smin <- 2
  df$eor <- 90
  df$x.err <- 5
  df$y.err <- 5

  expect_error(addMeasurementError(df),
               "Cannot provide both error ellipse and x- and y-axis error terms")
})

test_that("addMeasurementError correctly handles a mix of valid KF and LS errors in the same dataset", {
  df <- get_base_data()

  # Initialize all error columns with NAs
  df$smaj <- NA_real_
  df$smin <- NA_real_
  df$eor  <- NA_real_
  df$x.err <- NA_real_
  df$y.err <- NA_real_

  # Assign KF errors to the first 5 rows
  df$smaj[1:5] <- 5
  df$smin[1:5] <- 2
  df$eor[1:5]  <- 1.5

  # Assign LS/GPS errors to the last 5 rows
  df$x.err[6:10] <- 3
  df$y.err[6:10] <- 3

  set.seed(789)
  res <- suppressWarnings(addMeasurementError(df))

  # 1. Check that ALL rows received measurement error (x != mu.x)
  expect_false(any(res$x == res$mu.x))
  expect_false(any(res$y == res$mu.y))

  # 2. Check that no NAs were introduced into the final coordinates
  expect_true(all(!is.na(res$x)))
  expect_true(all(!is.na(res$y)))

  # 3. Verify the original error columns were perfectly preserved
  # (KF rows should still have NA for LS, and vice versa)
  expect_equal(res$smaj[1:5], rep(5, 5))
  expect_true(all(is.na(res$smaj[6:10])))

  expect_equal(res$x.err[6:10], rep(3, 5))
  expect_true(all(is.na(res$x.err[1:5])))
})

# --- Tests: Class Handling (dataLangevin / simLangevin) ---

test_that("addMeasurementError handles dataLangevin objects from formatData", {
  raw_df <- data.frame(
    id = "A",
    date = as.POSIXct("2024-01-01 10:00:00", tz = "UTC") + (1:5) * 3600,
    lon = c(0, 1, 2, 3, 4),
    lat = c(0, 1, 2, 3, 4)
  )

  fmt_dat <- suppressWarnings(formatData(raw_df, coord = c("lon", "lat")))

  expect_s3_class(fmt_dat, "dataLangevin")

  # Important Note: The Rcpp backend inherently relies on 'mu.x' and 'mu.y' columns existing
  # to extract the mean of the distribution, regardless of the `coord` argument passed.
  # Therefore, when adding simulated error to a raw dataLangevin object, we must explicitly map them.
  fmt_dat$mu.x <- fmt_dat$x
  fmt_dat$mu.y <- fmt_dat$y

  set.seed(789)
  res <- suppressWarnings(addMeasurementError(fmt_dat, measurementError = list(x.sd = 2, y.sd = 2)))

  expect_false(isTRUE(all.equal(res$x, res$mu.x)))
  expect_s3_class(res, "dataLangevin")
})

test_that("addMeasurementError handles simLangevin objects from simLangevin.default", {
  # EXPANDED EXTENT to prevent out-of-bounds C++ errors
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-1000, 1000, -1000, 1000))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  p <- list(beta = c(0.5), sigma = 5, gamma = 0.5)

  # Simulate pure tracks with no observation error
  sim_dat <- suppressMessages(simLangevin(model = "underdamped", par = p, spatialCovs = spatialCovs,
                                          nbAnimals = 1, obsPerAnimal = 10, timeStep = 1,
                                          initialPosition = c(50, 50)))

  expect_s3_class(sim_dat, "simLangevin")

  # Inject observation error post-hoc
  set.seed(123)
  res <- suppressWarnings(addMeasurementError(sim_dat, measurementError = list(x.sd = 2, y.sd = 2)))

  # Verify true locations were untouched while observed locations diverged
  expect_equal(res$mu.x, sim_dat$mu.x)
  expect_false(isTRUE(all.equal(res$x, res$mu.x)))

  # Ensure the class attributes weren't stripped
  expect_s3_class(res, "simLangevin")
})

test_that("addMeasurementError handles simLangevin objects from simLangevin.fitLangevin", {

  # 1. Create Dummy Data and Raster
  raw_df <- data.frame(
    id = rep("A", 10),
    date = as.POSIXct("2024-01-01 10:00:00", tz = "UTC") + (1:10) * 3600,
    x = seq(40, 50, length.out = 10),
    y = seq(40, 50, length.out = 10)
  )
  fmt_dat <- suppressWarnings(formatData(raw_df))

  # MASSIVE EXTENT to accommodate the 36,000-step random walk
  # triggered by POSIXct numeric conversion in the simulation
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-100000, 100000, -100000, 100000))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2. Fit Dummy Model
  # We set sigma to a tiny value (0.01) so the animal barely moves during the simulation
  fit <- suppressMessages(suppressWarnings(
    fitLangevin(data = fmt_dat, spatialCovs = spatialCovs,
                par = list(sigma = 0.01, gamma = 0.5), silent = TRUE)
  ))

  # 3. Simulate Forward Unconditionally
  sim_fit_dat <- suppressWarnings(suppressMessages(
    simLangevin(fit, data = fmt_dat, spatialCovs = spatialCovs, timeStep = "30 mins", conditional = FALSE)
  ))

  expect_s3_class(sim_fit_dat, "simLangevin")

  # 4. Generate New Post-hoc Errors (Kalman Filter this time)
  set.seed(456)
  res <- suppressWarnings(addMeasurementError(sim_fit_dat, measurementError = list(smaj.sd = 2, smin.sd = 1, eor.lim = c(0, 180))))

  # Verify divergence
  expect_false(isTRUE(all.equal(res$x, res$mu.x)))
  # Verify error components were successfully calculated and populated
  expect_true(all(!is.na(res$smaj)))
  expect_s3_class(res, "simLangevin")
})

# --- Tests: EMF Dataframe Errors ---

test_that("addMeasurementError catches invalid EMF dataframe structures", {
  df <- get_base_data()

  # Missing the 'emf.y' column
  bad_emf_cols <- data.frame(
    lc = c("3", "2"),
    emf.x = c(1, 2),
    prob = c(0.5, 0.5)
  )

  # Probabilities sum to 1.1 instead of 1.0
  bad_emf_prob <- data.frame(
    lc = c("3", "2"),
    emf.x = c(1, 2),
    emf.y = c(1, 2),
    prob = c(0.5, 0.6)
  )

  expect_error(
    addMeasurementError(df, measurementError = bad_emf_cols),
    "must contain columns 'lc', 'emf.x', 'emf.y', and 'prob'"
  )

  expect_error(
    addMeasurementError(df, measurementError = bad_emf_prob),
    "must sum to 1"
  )
})

test_that("addMeasurementError correctly applies EMF dataframe probabilities", {
  df <- get_base_data()

  # A valid custom EMF dataframe
  valid_emf <- data.frame(
    lc = c("3", "2", "1"),
    emf.x = c(1.0, 1.5, 3.0),
    emf.y = c(1.0, 1.5, 3.0),
    prob = c(0.5, 0.3, 0.2)
  )

  set.seed(42)
  res <- suppressWarnings(addMeasurementError(df, measurementError = valid_emf))

  # Verify the new columns were created
  expect_true("lc" %in% names(res))
  expect_true("x.err" %in% names(res))
  expect_true("y.err" %in% names(res))

  # Verify that Argos KF error columns are cleanly set to NA
  expect_true(all(is.na(res$smaj)))
  expect_true(all(is.na(res$smin)))
  expect_true(all(is.na(res$eor)))

  # Verify that the location classes drawn are strictly from our provided dataframe
  expect_true(all(res$lc %in% c("3", "2", "1")))

  # Verify errors were actually generated (no NAs in the LS/GPS columns)
  expect_false(any(is.na(res$x.err)))
  expect_false(any(is.na(res$y.err)))

  # Verify coordinates shifted
  expect_false(isTRUE(all.equal(res$x, res$mu.x)))
  expect_false(isTRUE(all.equal(res$y, res$mu.y)))
})
