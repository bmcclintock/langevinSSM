# tests/testthat/test_fitLangevin.R

# --- Helpers ---
get_valid_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  return(r)
}

get_valid_dataLangevin <- function() {
  df <- data.frame(
    id = as.factor(rep("A", 5)), # FIXED: Must be a factor to mimic formatData
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 11:00:00", "2023-01-01 12:00:00",
                        "2023-01-01 13:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    dt = rep(1, 5),
    x = c(40, 45, 50, 55, 60),
    y = c(40, 45, 50, 55, 60),
    x.sd = rep(1, 5),
    y.sd = rep(1, 5),
    smaj = rep(NA, 5),
    smin = rep(NA, 5),
    eor = rep(NA, 5),
    lc = as.factor(rep("G", 5)) # FIXED: Also made factor for consistency
  )
  attr(df, "time.unit") <- "hours"
  class(df) <- append("dataLangevin", class(df))
  return(df)
}

get_valid_par <- function() {
  list(beta = c(0.5), sigma = 5, gamma = 0.5)
}

test_that("fitLangevin rejects unformatted data", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # Just a regular data frame, not run through formatData()
  bad_data <- data.frame(id = "A", x = 50, y = 50)

  expect_error(fitLangevin(data = bad_data, spatialCovs = r, par = p),
               "is not formatted as a 'dataLangevin' object")
})

test_that("fitLangevin catches missing coordinate columns", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Tell fitLangevin to look for lon/lat instead of x/y
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, coord = c("lon", "lat")),
               "coord not found in data")
})

test_that("fitLangevin enforces smoothGradient rules", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Invalid npoints
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, npoints = 5),
               "npoints must be 4 or 8")

  # Invalid curweight (must be >= 0 and < 1)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, curweight = 1.5),
               "curweight must be >=0 and <1")
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, curweight = -0.1),
               "curweight must be >=0 and <1")

  # Invalid zetaScale
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, zetaScale = -2),
               "zetaScale must be >0")
})

test_that("fitLangevin catches spatial overlap errors from prepareRaster", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Push coordinates completely outside the 0-100 raster extent
  dat$x <- c(200, 205, 210, 215, 220)

  # This tests that fitLangevin correctly passes the data to prepareRaster
  # and fast-fails BEFORE attempting any slow TMB evaluations
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p),
               "overlap with 'spatialCovs'")
})

test_that("fitLangevin validates parameter lengths and map structures", {
  r <- list(habitat = get_valid_raster()) # 1 covariate
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # beta length mismatch (provided 2 coefficients for 1 covariate)
  p_bad_beta <- p
  p_bad_beta$beta <- c(0.5, -0.2)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_beta),
               "beta")

  # tau length mismatch (provided a single scalar instead of a 2-vector)
  p_bad_tau <- p
  p_bad_tau$tau <- 1.5
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_tau),
               "tau")

  # rho_o bounds check (correlation must theoretically be bounded, though TMB maps it)
  p_bad_rho <- p
  p_bad_rho$rho_o <- 1.5
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_rho),
               "rho_o")

  nonsense_p <- list(beta=0,gamma=0.5,sigma=1, nonsense=123)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = nonsense_p),
               "par")

  # map mismatch
  # If a user tries to map out tau (which is intrinsically length 2) but only provides a length-1 factor
  bad_map <- list(tau = as.factor(1))
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, map = bad_map),
               "map")

  # If a user tries to map out nonsense
  nonsense_map <- list(nonsense = as.factor(1))
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, map = nonsense_map),
               "map")
})

test_that("fitLangevin rejects observation parameters when data lacks corresponding errors", {
  r <- list(habitat = get_valid_raster())
  dat <- get_valid_dataLangevin() # This mock has x.sd/y.sd, but NO smaj/smin/eor
  p <- get_valid_par()

  # attempt to pass psi (should fail because there is no ellipse data)
  p_bad_psi <- p
  p_bad_psi$psi <- 1.5
  expect_error(suppressMessages(fitLangevin(data = dat, spatialCovs = r, par = p_bad_psi)),
               "error ellipse observations")

  # map psi should also fail
  expect_error(suppressMessages(fitLangevin(data = dat, spatialCovs = r, par = p, map = list(psi = as.factor(NA)))),
               "error ellipse observations")

  # now strip LS/GPS data to test tau/rho rejection
  dat_no_err <- dat
  dat_no_err$x.sd <- NA
  dat_no_err$y.sd <- NA

  p_bad_tau <- p
  p_bad_tau$tau <- c(1, 1)
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p_bad_tau)),
               "standard deviation observations")

  # map rho_o should fail
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p, map = list(rho_o = as.factor(NA)))),
               "standard deviation observations")
})

test_that("mapDuplicatedTimes correctly maps duplicate states and preserves user maps", {
  # create mock internal data structure
  dat <- list(
    ID = c("A", "A", "A", "B", "B"),
    Y = matrix(NA, nrow=2, ncol=5),
    dt = c(0, 1, 0, 0, 1) # Note: dt[3] is a duplicate for A
  )

  # mock parameters (5 locations = 10 elements)
  par <- list(
    mu = matrix(0, nrow=2, ncol=5),
    vel = matrix(0, nrow=2, ncol=5)
  )
  re <- c("mu", "vel")

  # mock a user map where the FIRST location (elements 1 and 2) is fixed to NA
  user_map <- list(
    mu = factor(c(NA, NA, 3:10))
  )

  # run the mapping function
  new_map <- langevinSSM:::mapDuplicatedTimes(dat, map = user_map, par = par, re = re)

  # --- Assertions ---
  mu_map_chr <- as.character(new_map$mu)

  # user map preservation
  # The first location should still be NA
  expect_true(is.na(mu_map_chr[1]) && is.na(mu_map_chr[2]))

  # duplicated Time Mapping (A)
  # dat$dt[3] == 0, so location 3 (elements 5,6) should map to location 2 (elements 3,4)
  expect_equal(mu_map_chr[5], mu_map_chr[3])
  expect_equal(mu_map_chr[6], mu_map_chr[4])

  # track boundary protection (B)
  # dat$dt[4] == 0, but this is the FIRST observation of track B.
  # It should NOT map to the last observation of track A (elements 5,6).
  expect_true(mu_map_chr[7] != mu_map_chr[5])

  # velocity Mapping
  # Since vel wasn't mapped by the user, it should have been initialized and mapped
  expect_false(is.null(new_map$vel))
  vel_map_chr <- as.character(new_map$vel)
  expect_equal(vel_map_chr[5], vel_map_chr[3]) # duplicate mapping applied to vel
})

test_that("mapDuplicatedTimes generates strict sequential TMB factors (lexicographical check)", {
  # 12 locations = 24 parameters (pushes past the "9" threshold to test alphabetical sorting)
  dat <- list(
    ID = rep("A", 12),
    Y = matrix(NA, nrow=2, ncol=12),
    dt = c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1) # duplicates at indices 6 and 11
  )
  par <- list(
    mu = matrix(0, nrow=2, ncol=12),
    vel = matrix(0, nrow=2, ncol=12)
  )
  re <- c("mu", "vel")

  # Create a map with NAs to verify NA preservation and sequential numbering
  user_map <- list(
    mu = factor(c(NA, NA, NA, NA, 5:24)) # Locations 1 and 2 are fixed
  )

  new_map <- langevinSSM:::mapDuplicatedTimes(dat, map = user_map, par = par, re = re)

  # TMB evaluates maps based on their internal integer representation
  mu_int <- as.integer(new_map$mu)

  # 1. Check NAs are preserved exactly where they belong
  expect_true(all(is.na(mu_int[1:4])))

  # 2. Extract the non-NA integers in the exact order they appear
  valid_mu <- mu_int[!is.na(mu_int)]

  # 3. Get the unique elements in the order of their first appearance
  first_appearances <- unique(valid_mu)

  # 4. TMB requires these to be strictly sequential (1, 2, 3...).
  # If R sorted the factor alphabetically, "10" would become an earlier factor level than "5",
  # causing this expectation to fail.
  expect_equal(first_appearances, 1:length(first_appearances))

  # 5. Check specific duplicate mapping
  # dt[6] == 0, so location 6 (elements 11, 12) maps to location 5 (elements 9, 10)
  expect_equal(mu_int[11], mu_int[9])
  expect_equal(mu_int[12], mu_int[10])
})

test_that("fitLangevin successfully fits end-to-end with duplicated times", {

  # 1. Simulate Dummy Spatial Covariate
  r <- terra::rast(nrows = 20, ncols = 20, ext = c(0, 200, 0, 200))
  terra::values(r) <- runif(terra::ncell(r)) # random habitat
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2. Simulate Track Data
  n_obs <- 30
  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dates <- start_time + (1:n_obs) * 3600 # 1 observation per hour

  set.seed(123)
  x <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)
  y <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)

  dat <- data.frame(
    id = rep("A", n_obs),
    date = dates,
    x = x,
    y = y,
    lc = rep("G", n_obs),
    x.sd = rep(3, n_obs), # 3m GPS error
    y.sd = rep(3, n_obs)
  )

  # 3. Inject Edge Cases
  # A. Hard Duplicate: Two actual observations logged at the exact same time
  dup_row <- dat[15, ]
  dup_row$x <- dup_row$x + 1.5
  dup_row$y <- dup_row$y - 0.5
  dat <- rbind(dat, dup_row)
  dat <- dat[order(dat$date), ] # Sort chronological

  # B. Soft Duplicate: Request a prediction exactly on top of an existing obs
  # Because we added a row at 15, the original 25th observation is now row 26
  pt <- data.frame(id = "A", date = dates[25])

  # 4. Format Data (Expect our custom warning)
  expect_warning(
    fmt_dat <- formatData(dat, predTimes = pt, time.unit = "hours"),
    "Duplicated times with observed"
  )

  expect_true(any(fmt_dat$dt == 0))

  # 5. Fit the Model
  init_par <- list(beta = 0, sigma = 5, gamma = 1)

  # We suppress warnings here to ignore the "NaNs produced" Hessian warning
  # that inherently occurs when fitting perfectly random noise, as well as
  # the 'dt < 1e-6' warning we explicitly programmed in.
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "underdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      silent = TRUE
    )
  }))

  # 6. Assertions
  expect_s3_class(fit, "fitLangevin")

  # Should converge (0) or hit max iterations safely (1), but not crash
  expect_true(fit$convergence %in% c(0, 1))

  # Verify the hard duplicate mapped correctly (rows 15 and 16)
  mu_x_15 <- fit$estimates$random$mu$est[15, "mu.x"]
  mu_x_16 <- fit$estimates$random$mu$est[16, "mu.x"]
  expect_equal(mu_x_15, mu_x_16)

  # Verify the SE extraction logic correctly expanded the SEs to both rows
  se_x_15 <- fit$estimates$random$mu$se[15, "mu.x"]
  se_x_16 <- fit$estimates$random$mu$se[16, "mu.x"]
  expect_equal(se_x_15, se_x_16)

  # Verify the soft duplicate (predTimes) was safely handled
  # Total rows should be 30 original + 1 hard duplicate = 31 (predTime collapsed)
  expect_equal(nrow(fit$estimates$random$mu$est), 31)
})

test_that("fitLangevin respects user-defined map with duplicated times", {

  # 1. Simulate Dummy Spatial Covariate
  r <- terra::rast(nrows = 20, ncols = 20, ext = c(0, 200, 0, 200))
  terra::values(r) <- runif(terra::ncell(r))
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2. Simulate Track Data (30 obs)
  n_obs <- 30
  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dates <- start_time + (1:n_obs) * 3600

  set.seed(456)
  x <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)
  y <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)

  dat <- data.frame(
    id = rep("A", n_obs),
    date = dates,
    x = x,
    y = y,
    lc = rep("G", n_obs),
    x.sd = rep(3, n_obs),
    y.sd = rep(3, n_obs)
  )

  # 3. Inject a Hard Duplicate (at index 10)
  dup_row <- dat[10, ]
  dup_row$x <- dup_row$x + 1.0
  dup_row$y <- dup_row$y - 1.0
  dat <- rbind(dat, dup_row)
  dat <- dat[order(dat$date), ] # Total 31 rows

  # 4. Format Data
  expect_warning(
    fmt_dat <- formatData(dat, time.unit = "hours"),
    "Duplicated times with observed"
  )

  # 5. Build Initial Parameters and Map
  N <- nrow(fmt_dat) # 31 observations

  # initialValues expects N rows and 2 columns for user inputs
  init_mu <- as.matrix(fmt_dat[, c("x", "y")])
  init_vel <- matrix(0, nrow = N, ncol = 2)

  # We want to freeze the 1st, 15th (middle), and 31st (last) locations to explicit values
  fixed_obs_indices <- c(1, 15, 31)

  # Inject obvious, weird values so we can easily verify they didn't move
  init_mu[fixed_obs_indices, 1] <- c(999, 888, 777)  # X values
  init_mu[fixed_obs_indices, 2] <- c(111, 222, 333)  # Y values

  init_vel[fixed_obs_indices, 1] <- c(10, 20, 30)    # X velocity
  init_vel[fixed_obs_indices, 2] <- c(-10, -20, -30) # Y velocity

  init_par <- list(
    beta = 0, sigma = 5, gamma = 1,
    mu = init_mu, vel = init_vel
  )

  # Create the user map to fix these indices
  # TMB internally evaluates parameters as a 1D column-major vector of its 2xN format:
  # (x1, y1, x2, y2 ...). So the fixed elements correspond to:
  # Obs 1:  indices 1 & 2
  # Obs 15: indices 29 & 30
  # Obs 31: indices 61 & 62
  mu_map <- 1:(2 * N)
  fix_elements <- c(1, 2, 29, 30, 61, 62)
  mu_map[fix_elements] <- NA

  user_map <- list(
    mu = factor(mu_map),
    vel = factor(mu_map)
  )

  # 6. Fit the Model
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "underdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      map = user_map,
      silent = TRUE
    )
  }))

  # 7. Assertions
  expect_s3_class(fit, "fitLangevin")
  expect_true(fit$convergence %in% c(0, 1))

  # A. Check Duplicate Logic still held
  expect_equal(
    fit$estimates$random$mu$est[10, "mu.x"],
    fit$estimates$random$mu$est[11, "mu.x"]
  )

  # B. Check Fixed Values remained exactly at their initialized state
  # Obs 1
  expect_equal(fit$estimates$random$mu$est[1, "mu.x"], 999)
  expect_equal(fit$estimates$random$mu$est[1, "mu.y"], 111)
  expect_equal(fit$estimates$random$vel$est[1, "vel.x"], 10)

  # Obs 15
  expect_equal(fit$estimates$random$mu$est[15, "mu.x"], 888)
  expect_equal(fit$estimates$random$mu$est[15, "mu.y"], 222)
  expect_equal(fit$estimates$random$vel$est[15, "vel.x"], 20)

  # Obs 31
  expect_equal(fit$estimates$random$mu$est[31, "mu.x"], 777)
  expect_equal(fit$estimates$random$mu$est[31, "mu.y"], 333)
  expect_equal(fit$estimates$random$vel$est[31, "vel.x"], 30)

  # C. Check Standard Errors for fixed values
  fixed_ses <- c(
    fit$estimates$random$mu$se[1, "mu.x"],
    fit$estimates$random$mu$se[15, "mu.x"],
    fit$estimates$random$mu$se[31, "mu.x"]
  )
  expect_true(all(is.na(fixed_ses)))
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
  dat$x.sd <- 2
  dat$y.sd <- 2
  dat$lc <- "G"
  dat$id <- "A"

  # 3. Inject Edge Cases

  # A. Hard Duplicate 1 (With observation error)
  dup_row1 <- dat[5, ]
  dup_row1$x <- dup_row1$x + 1.0 # Slight GPS jitter
  dat <- rbind(dat, dup_row1)

  # B. Hard Duplicate 2 (Known location, NO observation error)
  # We set x.sd and y.sd to NA to trigger the "perfect track" logic.
  dat$x.sd[20] <- NA
  dat$y.sd[20] <- NA
  dup_row2 <- dat[20, ]
  dat <- rbind(dat, dup_row2)

  dat <- dat[order(dat$date), ]

  # C. True Interpolation Point (NA coordinates)
  pt <- data.frame(id = dat$id[1], date = dat$date[10] + 1800)

  fmt_dat <- suppressWarnings(formatData(dat, predTimes = pt, time.unit = "hours"))

  # 4. Fit Model with calcOSA = TRUE
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "underdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      calcOSA = TRUE,
      silent = TRUE
    )
  }))

  # 5. Assertions
  expect_s3_class(fit, "fitLangevin")
  expect_true(is.data.frame(fit$osa))

  # Total rows should match formatted data exactly
  expect_equal(nrow(fit$osa), nrow(fmt_dat))

  # 1. The interpolation point (NA coordinates) should have NA residuals
  na_obs_idx <- which(is.na(fmt_dat$x))
  expect_true(length(na_obs_idx) == 1)
  expect_true(all(is.na(fit$osa$residual.x[na_obs_idx])))

  # 2. The first observation of EVERY track should have NA residuals
  first_obs_idx <- which(!duplicated(fmt_dat$id))
  expect_true(all(is.na(fit$osa$residual.x[first_obs_idx])))

  # 3. Known Locations (No observation error)
  # Because these have no error, TMB skips their observation likelihood (isd=0).
  # With no likelihood to evaluate, their residuals MUST be NA.
  known_obs_idx <- which(is.na(fmt_dat$x.sd) & !is.na(fmt_dat$x))
  expect_true(length(known_obs_idx) == 2) # The original + the duplicate
  expect_true(all(is.na(fit$osa$residual.x[known_obs_idx])))

  # 4. ALL OTHER real observations
  # This includes the hard duplicate WITH observation error! Since it has a valid
  # observation variance, TMB successfully computes a numeric residual for it.
  real_obs_idx <- which(!is.na(fmt_dat$x))

  # Filter out the first observation and the known (error-free) observations
  obs_to_check <- setdiff(real_obs_idx, c(first_obs_idx, known_obs_idx))

  expect_true(all(!is.na(fit$osa$residual.x[obs_to_check])))
  expect_true(all(!is.na(fit$osa$residual.y[obs_to_check])))
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
  dat$x.sd <- 2
  dat$y.sd <- 2
  dat$lc <- "G"
  dat$id <- "A"

  # 3. Inject Edge Cases

  # A. Hard Duplicate 1 (With observation error)
  dup_row1 <- dat[5, ]
  dup_row1$x <- dup_row1$x + 1.0 # Slight GPS jitter
  dat <- rbind(dat, dup_row1)

  # B. Hard Duplicate 2 (Known location, NO observation error)
  # We set x.sd and y.sd to NA to trigger the "perfect track" logic.
  dat$x.sd[20] <- NA
  dat$y.sd[20] <- NA
  dup_row2 <- dat[20, ]
  dat <- rbind(dat, dup_row2)

  dat <- dat[order(dat$date), ]

  # C. True Interpolation Point (NA coordinates)
  pt <- data.frame(id = dat$id[1], date = dat$date[10] + 1800)

  fmt_dat <- suppressWarnings(formatData(dat, predTimes = pt, time.unit = "hours"))

  # 4. Fit Model with calcOSA = TRUE
  suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "overdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      calcOSA = TRUE,
      silent = TRUE
    )
  }))

  # 5. Assertions
  expect_s3_class(fit, "fitLangevin")
  expect_true(is.data.frame(fit$osa))

  # Total rows should match formatted data exactly
  expect_equal(nrow(fit$osa), nrow(fmt_dat))

  # 1. The interpolation point (NA coordinates) should have NA residuals
  na_obs_idx <- which(is.na(fmt_dat$x))
  expect_true(length(na_obs_idx) == 1)
  expect_true(all(is.na(fit$osa$residual.x[na_obs_idx])))

  # 2. Known Locations (No observation error) MUST be NA
  known_obs_idx <- which(is.na(fmt_dat$x.sd) & !is.na(fmt_dat$x))
  expect_true(length(known_obs_idx) == 2) # The original + the duplicate
  expect_true(all(is.na(fit$osa$residual.x[known_obs_idx])))

  # 3. The first observation of EVERY track MUST be NA
  # (Because we explicitly conditioned on it in fitLangevin.R!)
  first_obs_idx <- which(!duplicated(fmt_dat$id))
  expect_true(all(is.na(fit$osa$residual.x[first_obs_idx])))
  expect_true(all(is.na(fit$osa$residual.y[first_obs_idx])))

  # 4. ALL OTHER real observations (including the hard duplicate WITH error)
  real_obs_idx <- which(!is.na(fmt_dat$x))

  # Filter out the known (error-free) observations and the first observation
  obs_to_check <- setdiff(real_obs_idx, c(first_obs_idx, known_obs_idx))

  expect_true(all(!is.na(fit$osa$residual.x[obs_to_check])))
  expect_true(all(!is.na(fit$osa$residual.y[obs_to_check])))
})
