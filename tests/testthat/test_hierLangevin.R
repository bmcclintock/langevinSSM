# tests/testthat/test_hierLangevin.R

# ==============================================================================
# 1. MOCK DATA GENERATION HELPER
# ==============================================================================
# Generates a lightweight, valid fitLangevin object to bypass expensive Stage I SDE
# optimization during unit tests. Includes the barrier and scaleFactor arguments.
make_mock_fit <- function(id, param_names, mu_vals, sd_val = 0.1, barrier = NULL, lambda = NULL, scaleFactor = 1) {
  p <- length(param_names)

  # Create a positive-definite covariance matrix
  cov_mat <- diag(sd_val^2, p)
  dimnames(cov_mat) <- list(param_names, param_names)

  structure(
    list(
      estimates = list(
        working = data.frame(Estimate = mu_vals, row.names = param_names)
      ),
      covariance = list(
        working = cov_mat
      ),
      conditions = list(
        model = "underdamped",
        ids = id,
        barrier = barrier,
        lambda = lambda,
        scaleFactor = scaleFactor
      )
    ),
    class = "fitLangevin"
  )
}

# Setup a valid list of 3 mock individual fits, explicitly including barrier/sf metadata
p_names <- c("beta_cov1", "log_sigma", "log_gamma")
valid_fit_list <- list(
  ind_1 = make_mock_fit("ind_1", p_names, c(-1.0, 1.5, -0.5), barrier = "coast", lambda = 100, scaleFactor = 1000),
  ind_2 = make_mock_fit("ind_2", p_names, c(-1.2, 1.6, -0.4), barrier = "coast", lambda = 100, scaleFactor = 1000),
  ind_3 = make_mock_fit("ind_3", p_names, c(-0.8, 1.4, -0.6), barrier = "coast", lambda = 100, scaleFactor = 1000)
)


# ==============================================================================
# 2. UNIT TESTS
# ==============================================================================

test_that("hierLangevin input validation and filtering works correctly", {

  # Fails if not a list, or if the user passed a single fitLangevin object directly
  expect_error(hierLangevin(valid_fit_list[[1]]), "'fit_list' must be a list containing")

  # Fails if list is empty
  expect_error(hierLangevin(list()), "'fit_list' must be a list containing")

  # Fails if less than 2 valid individuals
  expect_error(hierLangevin(valid_fit_list[1]), "requires at least 2 successful individual fits")

  # Drops invalid list elements with a warning
  mixed_list <- valid_fit_list
  mixed_list$ind_4 <- "this is not a fitLangevin object"
  mixed_list$ind_5 <- structure(list(), class = "try-error") # Simulates a failed try() block

  warns <- capture_warnings(fit_dropped <- hierLangevin(mixed_list, silent = TRUE))
  expect_true(any(grepl("Dropped 2 element\\(s\\) from 'fit_list'", warns)))
  expect_s3_class(fit_dropped, "hierLangevin")
})


test_that("hierLangevin detects structural base mismatches between Stage I fits", {

  # Mismatch: Different parameter names
  bad_names_list <- valid_fit_list
  bad_names_list$ind_2 <- make_mock_fit("ind_2", c("beta_cov2", "log_sigma", "log_gamma"), c(-1,1,1), barrier = "coast", lambda = 100, scaleFactor = 1000)
  expect_error(
    hierLangevin(bad_names_list),
    "parameter names for individual 'ind_2' do not match"
  )

  # Mismatch: Missing Stage I covariance matrix (calcSE was FALSE)
  no_se_list <- valid_fit_list
  no_se_list$ind_3$covariance$working <- NULL
  expect_error(
    hierLangevin(no_se_list),
    "Working-scale covariance matrix unavailable for individual 'ind_3'"
  )
})


test_that("hierLangevin strictly enforces identical barrier, lambda, and scaleFactor arguments", {

  # Alter the barrier name
  bad_barrier <- valid_fit_list
  bad_barrier$ind_2$conditions$barrier <- "different_coast"
  expect_error(
    hierLangevin(bad_barrier),
    "Barrier argument for individual 'ind_2' does not match individual 'ind_1'"
  )

  # Alter the lambda penalty
  bad_lambda <- valid_fit_list
  bad_lambda$ind_2$conditions$lambda <- 50
  expect_error(
    hierLangevin(bad_lambda),
    "Lambda penalty for individual 'ind_2' does not match individual 'ind_1'"
  )

  # Alter the scaleFactor
  bad_sf <- valid_fit_list
  bad_sf$ind_3$conditions$scaleFactor <- 500
  expect_error(
    hierLangevin(bad_sf),
    "scaleFactor argument for individual 'ind_3' does not match individual 'ind_1'"
  )
})


test_that("hierLangevin Maximum Penalized Likelihood (MPL) argument parsing works", {

  # MPL < 1 throws a warning because it pushes variance to 0.
  # Because this deliberately breaks the model, nlminb throws multiple NA/NaN
  # and convergence warnings. We capture all warnings to prevent test clutter.
  warns <- capture_warnings(hierLangevin(valid_fit_list, mpl = 0.5, silent = TRUE))
  expect_true(any(grepl("alpha < 1 applies a penalty", warns)))

  # MPL length mismatch throws an error
  expect_error(
    hierLangevin(valid_fit_list, mpl = c(2, 2)), # p = 3, so length 2 is invalid
    "must be a scalar or a vector of length 3"
  )

  # Valid MPL scalar passes silently
  expect_s3_class(hierLangevin(valid_fit_list, mpl = 1.5, silent = TRUE), "hierLangevin")

  # Valid MPL vector passes silently
  expect_s3_class(hierLangevin(valid_fit_list, mpl = c(2, 1.5, 2), silent = TRUE), "hierLangevin")
})


test_that("hierLangevin successfully optimizes and transfers metadata", {

  fit <- hierLangevin(valid_fit_list, silent = TRUE)

  # Check base properties
  expect_s3_class(fit, "hierLangevin")
  expect_true(fit$convergence == 0)
  expect_type(fit$objective, "double")

  # Verify metadata was properly grabbed from ind_1 and stored for downstream functions (e.g. getUD)
  expect_equal(fit$conditions$barrier, "coast")
  expect_equal(fit$conditions$lambda, 100)
  expect_equal(fit$conditions$scaleFactor, 1000)

  # Check Working Scale Estimates (Should have 2*p rows: means and sds)
  p <- length(p_names)
  expect_equal(nrow(fit$estimates$working), 2 * p)
  expect_true(all(c("Estimate", "Std. Error") %in% colnames(fit$estimates$working)))

  # Verify working scale rownames conform to mu_ and sd_ prefixes
  expected_working_rows <- c(paste0("mu_", p_names), paste0("sd_", p_names))
  expect_equal(rownames(fit$estimates$working), expected_working_rows)

  # Check Natural Scale Estimates (Should have p rows, and drop the 'log_' prefix)
  expect_equal(nrow(fit$estimates$natural), p)
  expected_nat_rows <- c("mu_beta_cov1", "mu_sigma", "mu_gamma")
  expect_equal(rownames(fit$estimates$natural), expected_nat_rows)

  # Check Random Effects / BLUPs
  expect_type(fit$estimates$random, "list")
  expect_equal(names(fit$estimates$random), c("beta_cov1", "sigma", "gamma"))

  # Check internal structure of a specific BLUP
  gamma_blup <- fit$estimates$random$gamma
  expect_s3_class(gamma_blup$est, "data.frame")
  expect_equal(nrow(gamma_blup$est), 3) # 3 individuals
  expect_equal(colnames(gamma_blup$est), c("id", "est"))
  expect_equal(as.character(gamma_blup$est$id), c("ind_1", "ind_2", "ind_3"))
})
