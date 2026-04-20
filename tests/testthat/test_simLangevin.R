# tests/testthat/test_simLangevin.R

# --- Helpers ---
get_valid_raster <- function() {
  # Expanded extent to prevent out-of-bounds errors during random walks
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-1000, 1000, -1000, 1000))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  return(r)
}

get_dynamic_raster <- function(type = "numeric") {
  # Expanded extent to prevent out-of-bounds errors during random walks
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-1000, 1000, -1000, 1000))
  terra::values(r) <- 1:100
  r_multi <- c(r, r)
  if (type == "numeric") {
    terra::time(r_multi) <- c(0, 10)
  } else {
    terra::time(r_multi) <- as.POSIXct(c("2023-01-01", "2023-01-02"), tz="UTC")
  }
  names(r_multi) <- c("habitat", "habitat")
  return(r_multi)
}

get_valid_par <- function() {
  list(beta = c(0.5), sigma = 5, gamma = 0.5)
}

test_that("simLangevin catches invalid numerical arguments", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = -5), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = Inf), "should be at least 1")

  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = -10), "should be at least 1")

  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = 0), "greater than zero")
  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = -0.5), "greater than zero")
})

test_that("simLangevin catches invalid model types", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  expect_error(simLangevin(model = "invalid_model", par = p, spatialCovs = r), "should be one of")
})

test_that("simLangevin returns a correctly structured dataLangevin object", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  res <- suppressMessages(simLangevin(model = "underdamped", par = p, spatialCovs = r,
                                      nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1,
                                      initialPosition = c(50, 50)))

  expect_s3_class(res, "dataLangevin")
  expect_true(all(c("id", "date", "dt", "x", "y", "mu.x", "mu.y", "vel.x", "vel.y") %in% names(res)))
  expect_equal(nrow(res), 5)

  # Remove gamma for the overdamped model!
  p_overdamped <- list(beta = p$beta, sigma = p$sigma)

  res_over <- suppressMessages(simLangevin(model = "overdamped", par = p_overdamped, spatialCovs = r,
                                           nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1,
                                           initialPosition = c(50, 50)))
  expect_true(all(c("mu.x", "mu.y") %in% names(res_over)))
  expect_false("vel.x" %in% names(res_over))
})

test_that("simLangevin validates natural scale parameters in par", {
  r <- list(habitat = get_valid_raster())

  expect_error(simLangevin(par = list(sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  expect_error(simLangevin(par = list(beta = c(0.5, -0.2), sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  expect_error(simLangevin(par = list(beta = c(0.5), gamma = 0.5), spatialCovs = r),
               "sigma")

  expect_error(simLangevin(model = "underdamped", par = list(beta = c(0.5), sigma = 5), spatialCovs = r),
               "gamma")

  expect_error(simLangevin(par = list(beta = c(0.5), sigma = 5, gamma = 0.5, tau = 1.5), spatialCovs = r),
               "tau")
})

test_that("simLangevin.default catches POSIXt times for dynamic covariates", {
  # dynamic raster with POSIXt times
  r <- list(habitat = get_dynamic_raster("POSIXt"))
  p <- get_valid_par()

  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = 5),
               "must be numeric \\(not POSIXt or Date\\)")
})

test_that("simLangevin.default enforces numeric temporal bounding limits", {
  # dynamic raster with numeric times (0 and 10)
  r <- list(habitat = get_dynamic_raster("numeric"))
  p <- get_valid_par()

  # Will simulate from 0 to (15-1)*1 = 14. Max time > 10, should fail.
  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = 15, timeStep = 1),
               "fall outside the temporal boundaries")

  # Will simulate from 0 to (5-1)*1 = 4. Max time < 10, should pass.
  res <- suppressMessages(simLangevin(model = "underdamped", par = p, spatialCovs = r,
                                      nbAnimals = 1, obsPerAnimal = 5, timeStep = 1,
                                      initialPosition = c(50, 50)))
  expect_s3_class(res, "dataLangevin")
})

test_that("simLangevin validates EMF dataframe structure and probabilities", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

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
    suppressMessages(simLangevin(par = p, spatialCovs = r, measurementError = bad_emf_cols)),
    "must contain columns 'lc', 'emf.x', 'emf.y', and 'prob'"
  )

  expect_error(
    suppressMessages(simLangevin(par = p, spatialCovs = r, measurementError = bad_emf_prob)),
    "must sum to 1"
  )
})

test_that("simLangevin correctly applies EMF dataframe probabilities to tracks", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # A valid custom EMF dataframe
  valid_emf <- data.frame(
    lc = c("3", "2", "1"),
    emf.x = c(1.0, 1.5, 3.0),
    emf.y = c(1.0, 1.5, 3.0),
    prob = c(0.5, 0.3, 0.2)
  )

  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  res <- suppressMessages(simLangevin(
    model = "underdamped",
    par = p,
    spatialCovs = r,
    nbAnimals = 1,
    obsPerAnimal = 50,
    timeStep = 0.1,
    initialPosition = c(50, 50),
    measurementError = valid_emf
  ))

  # Verify the new columns were preserved in the final output
  expect_true("lc" %in% names(res))
  expect_true("x.err" %in% names(res))
  expect_true("y.err" %in% names(res))

  # Verify that Argos KF error columns are cleanly set to NA
  expect_true(all(is.na(res$smaj)))
  expect_true(all(is.na(res$smin)))

  # Verify that the location classes drawn are strictly from our provided dataframe
  expect_true(all(res$lc %in% c("3", "2", "1")))

  # Verify errors were actually generated (no NAs in the LS/GPS columns)
  expect_false(any(is.na(res$x.err)))
  expect_false(any(is.na(res$y.err)))
})
