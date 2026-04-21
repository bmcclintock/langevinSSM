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
    id = as.factor(rep("A", 5)), # must be a factor to mimic formatData
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 11:00:00", "2023-01-01 12:00:00",
                        "2023-01-01 13:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    dt = rep(1, 5),
    x = c(40, 45, 50, 55, 60),
    y = c(40, 45, 50, 55, 60),
    x.err = rep(1, 5),
    y.err = rep(1, 5),
    smaj = rep(NA, 5),
    smin = rep(NA, 5),
    eor = rep(NA, 5),
    lc = as.factor(rep("G", 5)) # also made factor for consistency
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
  dat <- get_valid_dataLangevin() # This mock has x.err/y.err, but NO smaj/smin/eor
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
  dat_no_err$x.err <- NA
  dat_no_err$y.err <- NA

  p_bad_tau <- p
  p_bad_tau$tau <- c(1, 1)
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p_bad_tau)),
               "standard error observations")

  # map rho_o should fail
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p, map = list(rho_o = as.factor(NA)))),
               "standard error observations")
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

  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  x <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)
  y <- seq(50, 150, length.out = n_obs) + rnorm(n_obs, 0, 2)

  dat <- data.frame(
    id = rep("A", n_obs),
    date = dates,
    x = x,
    y = y,
    lc = rep("G", n_obs),
    x.err = rep(3, n_obs), # 3m GPS error
    y.err = rep(3, n_obs)
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

  class(dat) <- "data.frame"
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
    x.err = rep(3, n_obs),
    y.err = rep(3, n_obs)
  )

  # 3. Inject a Hard Duplicate (at index 10)
  dup_row <- dat[10, ]
  dup_row$x <- dup_row$x + 1.0
  dup_row$y <- dup_row$y - 1.0
  dat <- rbind(dat, dup_row)
  dat <- dat[order(dat$date), ] # Total 31 rows

  class(dat) <- "data.frame"

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

  # 6. Fit the Model (capture.output suppresses the boundsWarning cat() printout)
  capture.output(suppressMessages(suppressWarnings({
    fit <- fitLangevin(
      data = fmt_dat,
      model = "underdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      map = user_map,
      silent = TRUE
    )
  })))

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
  expect_true(all(fixed_ses == 0))
})

test_that("S3 methods for fitLangevin work", {

  # Build a lightweight, mock fitLangevin object with all necessary components
  mock_fit <- list(
    par = c(beta_hab = 1.5, sigma = 2),
    objective = 150.5,
    estimates = list(
      natural = data.frame(Estimate = c(1.5, 2.0), "Std. Error" = c(0.1, 0.2), row.names = c("beta_hab", "sigma"), check.names = FALSE),
      working = data.frame(Estimate = c(1.5, log(2)), "Std. Error" = c(0.1, 0.1), row.names = c("beta_hab", "log_sigma"), check.names = FALSE),
      random = list(
        mu = list(
          est = data.frame(id = c("A", "A"), date = as.POSIXct(c("2023-01-01", "2023-01-02"), tz="UTC"), mu.x = c(10, 20), mu.y = c(15, 25)),
          se = data.frame(id = c("A", "A"), date = as.POSIXct(c("2023-01-01", "2023-01-02"), tz="UTC"), mu.x = c(1, 1), mu.y = c(2, 2))
        )
      )
    ),
    covariance = list(
      natural = matrix(c(0.01, 0.005, 0.005, 0.04), nrow = 2, dimnames = list(c("beta_hab", "sigma"), c("beta_hab", "sigma")))
    ),
    signatures = list(data = list(nrow = 50)),
    convergence = 0,
    elapsedTime = c(user = 1, system = 0, elapsed = 1.5)
  )
  class(mock_fit) <- "fitLangevin"

  # 1. coef() and vcov()
  cf <- coef(mock_fit)
  expect_equal(length(cf), 2)
  expect_equal(cf["beta_hab"], c(beta_hab = 1.5))

  vc <- vcov(mock_fit)
  expect_true(is.matrix(vc))
  expect_equal(vc["sigma", "sigma"], 0.04)

  # 2. logLik(), AIC(), BIC()
  ll <- logLik(mock_fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -150.5)
  expect_equal(attr(ll, "df"), 2)
  expect_equal(attr(ll, "nobs"), 50)

  expect_equal(AIC(mock_fit), -2 * (-150.5) + 2 * 2)
  expect_equal(BIC(mock_fit), -2 * (-150.5) + log(50) * 2)

  # 3. confint() for Fixed Effects
  ci_fixed <- confint(mock_fit, type = "natural")
  expect_true(is.matrix(ci_fixed))
  expect_equal(colnames(ci_fixed), c("2.5 %", "97.5 %"))
  expect_true(ci_fixed["beta_hab", "2.5 %"] < 1.5)

  # Error catching in confint
  expect_error(confint(mock_fit, parm = "nonsense"), "Parameter\\(s\\) not found")

  # 4. confint() for Random Effects (wide layout, no point estimates)
  ci_mu <- confint(mock_fit, type = "mu")
  expect_true(is.data.frame(ci_mu))
  # Should have 6 columns: id, date, x bounds, y bounds
  expect_equal(colnames(ci_mu), c("id", "date", "mu.x_2.5%", "mu.x_97.5%", "mu.y_2.5%", "mu.y_97.5%"))
  expect_equal(nrow(ci_mu), 2)

  # Verify auto-correction of parm="mu" to type="mu"
  ci_mu_param_intercept <- confint(mock_fit, parm = "mu")
  expect_equal(ci_mu, ci_mu_param_intercept)

  # 5. print() output formatting
  out <- capture.output(print(mock_fit))
  expect_true(any(grepl("Habitat-Driven Langevin Diffusion Model", out)))
  expect_true(any(grepl("Parameter Estimates \\(Natural Scale\\):", out)))
  expect_true(any(grepl("beta_hab", out)))
})

test_that("summary, fitted, and numeric dates work correctly across S3 methods", {

  # Build mock fit with numeric dates and vel explicitly added
  mock_fit_num <- list(
    par = c(beta_hab = 1.5, sigma = 2),
    objective = 150.5,
    estimates = list(
      natural = data.frame(Estimate = c(1.5, 2.0), "Std. Error" = c(0.1, 0.2), row.names = c("beta_hab", "sigma"), check.names = FALSE),
      working = data.frame(Estimate = c(1.5, log(2)), "Std. Error" = c(0.1, 0.1), row.names = c("beta_hab", "log_sigma"), check.names = FALSE),
      random = list(
        mu = list(
          est = data.frame(id = c("A", "A"), date = c(10, 20), mu.x = c(10, 20), mu.y = c(15, 25)),
          se = data.frame(id = c("A", "A"), date = c(10, 20), mu.x = c(1, 1), mu.y = c(2, 2))
        ),
        vel = list(
          est = data.frame(id = c("A", "A"), date = c(10, 20), vel.x = c(0.1, 0.2), vel.y = c(0.3, 0.4)),
          se = data.frame(id = c("A", "A"), date = c(10, 20), vel.x = c(0.01, 0.01), vel.y = c(0.02, 0.02))
        )
      )
    ),
    covariance = list(
      natural = matrix(c(0.01, 0.005, 0.005, 0.04), nrow = 2, dimnames = list(c("beta_hab", "sigma"), c("beta_hab", "sigma")))
    ),
    signatures = list(data = list(nrow = 50)),
    convergence = 0,
    message = "relative convergence",
    elapsedTime = c(user = 1, system = 0, elapsed = 1.5)
  )
  class(mock_fit_num) <- "fitLangevin"

  # 1. fitted() method
  fit_mu <- fitted(mock_fit_num) # default is "mu"
  expect_true(is.data.frame(fit_mu))
  expect_true("date" %in% names(fit_mu))
  expect_true(is.numeric(fit_mu$date)) # Because we initialized with numeric
  expect_equal(fit_mu$date, c(10, 20))

  fit_vel <- fitted(mock_fit_num, type = "vel")
  expect_true("vel.x" %in% names(fit_vel))
  expect_true(is.numeric(fit_vel$date))

  # Error catching in fitted
  expect_error(fitted(mock_fit_num, type = "nonsense"), "should be one of")

  # 2. confint() with numeric dates
  ci_mu <- confint(mock_fit_num, type = "mu")
  expect_true("date" %in% names(ci_mu))
  expect_true(is.numeric(ci_mu$date))
  expect_equal(ci_mu$date, c(10, 20))

  ci_vel <- confint(mock_fit_num, type = "vel")
  expect_true("date" %in% names(ci_vel))
  expect_true(is.numeric(ci_vel$date))
  expect_equal(colnames(ci_vel), c("id", "date", "vel.x_2.5%", "vel.x_97.5%", "vel.y_2.5%", "vel.y_97.5%"))

  # 3. summary() method
  sum_fit <- summary(mock_fit_num)
  expect_s3_class(sum_fit, "summary.fitLangevin")

  # Ensure p-values were only calculated for habitat selection coefficients (beta)
  expect_true("Pr(>|z|)" %in% colnames(sum_fit$coef_beta))
  expect_false("Pr(>|z|)" %in% colnames(sum_fit$coef_process))
  expect_equal(rownames(sum_fit$coef_beta), "beta_hab")
  expect_equal(rownames(sum_fit$coef_process), "sigma")

  # 4. print.summary() method
  out_sum <- capture.output(print(sum_fit))
  expect_true(any(grepl("Habitat Selection Coefficients:", out_sum)))
  expect_true(any(grepl("Process & Observation Parameters (Natural Scale):", out_sum,fixed=TRUE)))
  expect_true(any(grepl("beta_hab", out_sum)))
  expect_true(any(grepl("sigma", out_sum)))
})

test_that("Out-of-bounds warnings trigger across fit and downstream methods", {
  # 1. Setup Data and Raster
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  dat <- data.frame(
    id = as.factor(c("A", "A")),
    date = as.POSIXct(c("2024-01-01 00:00:00", "2024-01-01 01:00:00"), tz="UTC"),
    dt = c(0, 1),
    x = c(5, NA),
    y = c(5, NA),
    lc = as.factor(c("G", "G")),
    x.err = c(1, NA),
    y.err = c(1, NA),
    smaj = as.numeric(c(NA, NA)),
    smin = as.numeric(c(NA, NA)),
    eor = as.numeric(c(NA, NA))
  )
  attr(dat, "time.unit") <- "hours"
  class(dat) <- c("dataLangevin", "data.frame")

  init_par <- list(
    beta = c(0), sigma = 1, gamma = 1,
    mu = matrix(c(5, 5, 9.5, 9.5), nrow = 2, byrow = TRUE),
    vel = matrix(c(0, 0, 0, 0), nrow = 2, byrow = TRUE)
  )

  user_map <- list(
    beta = factor(NA), sigma = factor(NA), gamma = factor(NA),
    mu = factor(c(1, NA, NA, NA)), vel = factor(rep(NA, 4))
  )

  # 2. Test Computation Functions (Expect Formal Warnings)
  # Chain expect_warning to absorb both the custom bounds warning
  # AND the TMB "empty summary" warning caused by the completely frozen mock map
  expect_warning(
    expect_warning(
      suppressMessages({
        fit <- fitLangevin(
          data = dat, model = "underdamped", spatialCovs = spatialCovs,
          par = init_par, map = user_map, silent = TRUE
        )
      }),
      "MODEL FIT LIKELY INVALID"
    ),
    "empty summary"
  )

  expect_true(fit$conditions$out_of_bounds)

  # residuals() throws the bounds warning FIRST, then an OSA warning due to the tiny mock dataset.
  expect_warning(
    expect_warning(
      suppressMessages(try(residuals(fit, dat, spatialCovs, run_tests = FALSE), silent=TRUE)),
      "MODEL FIT LIKELY INVALID"
    ),
    "OSA calculation failed"
  )

  expect_warning(
    suppressMessages(getUD(spatialCovs, fit = fit, log = TRUE, plot = FALSE)),
    "MODEL FIT LIKELY INVALID"
  )

  expect_warning(
    suppressMessages(plot(fit, spatialCovs = spatialCovs)),
    "MODEL FIT LIKELY INVALID"
  )

  # 3. Test Inspection Functions (Expect Printed Text via cat)
  out_print <- capture.output(print(fit))
  expect_true(any(grepl("WARNING: MODEL FIT LIKELY INVALID", out_print)))

  out_summary <- capture.output(print(summary(fit)))
  expect_true(any(grepl("WARNING: MODEL FIT LIKELY INVALID", out_summary)))
})
