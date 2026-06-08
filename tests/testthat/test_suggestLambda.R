# tests/testthat/test_suggestLambda.R

# ==============================================================================
# 1. SETUP MOCK DATA
# ==============================================================================

# Create a valid mock unconstrained overdamped fitLangevin object
mock_fit_od <- list(
  conditions = list(
    model = "overdamped",
    scaleFactor = 1,
    lambda = 0
  ),
  estimates = list(
    natural = data.frame(
      Estimate = c(2.0),
      row.names = c("sigma")
    )
  )
)
class(mock_fit_od) <- "fitLangevin"

# Create a valid mock unconstrained underdamped fitLangevin object
mock_fit_ud <- list(
  conditions = list(
    model = "underdamped",
    scaleFactor = 1,
    lambda = NULL # NULL is also valid for unconstrained
  ),
  estimates = list(
    natural = data.frame(
      Estimate = c(2.0, 0.5),
      row.names = c("sigma", "gamma")
    )
  )
)
class(mock_fit_ud) <- "fitLangevin"

# Create a mock CONSTRAINED fitLangevin object (should fail)
mock_fit_constrained <- mock_fit_ud
mock_fit_constrained$conditions$lambda <- 100

# ==============================================================================
# 2. TEST INPUT VALIDATION & USER ERRORS
# ==============================================================================

test_that("suggestLambda fails fast on invalid fit classes", {
  bad_fit <- unclass(mock_fit_od)

  expect_error(
    suggestLambda(bad_fit, max_dt = 3600),
    "'fit' must be provided and must be a fitLangevin object"
  )
})

test_that("suggestLambda fails if fit was constrained (lambda != 0)", {
  expect_error(
    suggestLambda(mock_fit_constrained, max_dt = 3600),
    "The provided 'fit' object must be an unconstrained model fitted with lambda = 0"
  )
})

test_that("suggestLambda fails fast on invalid max_dt values", {
  # Missing max_dt
  expect_error(
    suggestLambda(mock_fit_od),
    "'max_dt' must be provided as a single positive numeric value."
  )

  # Non-numeric
  expect_error(
    suggestLambda(mock_fit_od, max_dt = "1 hour"),
    "'max_dt' must be provided as a single positive numeric value."
  )

  # Negative value
  expect_error(
    suggestLambda(mock_fit_od, max_dt = -3600),
    "'max_dt' must be provided as a single positive numeric value."
  )

  # Vector of multiple values
  expect_error(
    suggestLambda(mock_fit_od, max_dt = c(1800, 3600)),
    "'max_dt' must be provided as a single positive numeric value."
  )
})

# ==============================================================================
# 3. TEST MATHEMATICAL CALCULATIONS
# ==============================================================================

test_that("suggestLambda accurately calculates lambda for overdamped models", {
  # For sigma = 2, max_dt = 3600:
  # lambda = 2 / (sigma^2 * max_dt) = 2 / (4 * 3600) = 2 / 14400 = 0.0001388889

  expected_lambda <- 2 / (4 * 3600)

  expect_message(
    calc_val <- suggestLambda(mock_fit_od, max_dt = 3600),
    "Suggested maximum barrier penalty"
  )

  expect_equal(calc_val, expected_lambda, tolerance = 1e-6)
})

test_that("suggestLambda accurately calculates lambda for underdamped models", {
  # For sigma = 2, gamma = 0.5, max_dt = 3600
  # num = (0.5^2) * (1 - exp(-0.5 * 3600)) = 0.25 * (1 - 0) = 0.25
  # den = (2^2) * (1 - exp(-1800) - (1800 * exp(-1800))) = 4 * (1 - 0 - 0) = 4
  # lambda = 0.25 / 4 = 0.0625

  expected_lambda <- 0.0625

  expect_message(
    calc_val <- suggestLambda(mock_fit_ud, max_dt = 3600),
    "Suggested maximum barrier penalty"
  )

  expect_equal(calc_val, expected_lambda, tolerance = 1e-6)
})
