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


# ==============================================================================
# 3. BARRIER STRENGTH (lambda) MODES
# ==============================================================================

make_lambda_fit <- function(id, lam_nat, sf = 1, mode = "estimate", se = 0.1, beta = -1, lse = 0.3) {
  if (mode == "estimate") {
    pn <- c("beta_cov1", "log_sigma", "log_gamma", "log_lambda")
    # working scale (coordinates / scaleFactor): log_sigma - log(sf), log_lambda + 2 log(sf)
    f <- make_mock_fit(id, pn, c(beta, 1.5 - log(sf), -0.5, log(lam_nat) + 2 * log(sf)), sd_val = se,
                       barrier = "coast", lambda = lam_nat, scaleFactor = sf)
    f$covariance$working["log_lambda", "log_lambda"] <- lse^2
  } else {
    f <- make_mock_fit(id, p_names, c(beta, 1.5 - log(sf), -0.5), sd_val = se, barrier = "coast",
                       lambda = if (mode == "tied") Inf else lam_nat, scaleFactor = sf)
  }
  f$conditions$lambda_mode <- mode
  f
}

test_that("hierLangevin pools log(lambda) on the natural scale when lambda = 'estimate'", {
  lam <- c(0.02, 0.05, 0.1, 0.2)
  b <- c(-1.0, -1.2, -0.8, -1.1)
  fl1 <- setNames(lapply(1:4, function(i) make_lambda_fit(paste0("ind_", i), lam[i], sf = 1, beta = b[i])), paste0("ind_", 1:4))
  fl10 <- setNames(lapply(1:4, function(i) make_lambda_fit(paste0("ind_", i), lam[i], sf = 10, beta = b[i])), paste0("ind_", 1:4))
  h1 <- hierLangevin(fl1, silent = TRUE)
  h10 <- hierLangevin(fl10, silent = TRUE)
  expect_s3_class(h1, "hierLangevin")
  expect_equal(h1$convergence, 0)
  expect_true(all(c("mu_log_lambda", "sd_log_lambda") %in% rownames(h1$estimates$working)))
  expect_true("mu_lambda" %in% rownames(h1$estimates$natural))
  expect_equal(names(h1$estimates$random), c("beta_cov1", "sigma", "gamma", "lambda"))
  expect_equal(h1$conditions$lambda, "estimate")
  expect_equal(h1$conditions$lambda_mode, "estimate")
  # pooled on the natural scale: mean of log(lambda) is near mean(log(lam)) and results do not depend on scaleFactor
  expect_equal(h1$estimates$working["mu_log_lambda", "Estimate"], mean(log(lam)), tolerance = 0.05)
  expect_equal(h1$estimates$working$Estimate, h10$estimates$working$Estimate, tolerance = 1e-6)
  expect_equal(h1$estimates$natural$Estimate, h10$estimates$natural$Estimate, tolerance = 1e-6)
  expect_equal(h1$estimates$natural["mu_lambda", "Estimate"], exp(h1$estimates$working["mu_log_lambda", "Estimate"]), tolerance = 1e-8)
  # individual lambdas shrink towards the population mean but keep their ordering
  ind_lam <- h1$estimates$random$lambda$est$est
  expect_true(all(ind_lam > 0))
  expect_equal(order(ind_lam), 1:4)
})

test_that("hierLangevin requires a common lambda mode, and a common fixed value", {
  fl <- list(a = make_lambda_fit("a", 0.05), b = make_lambda_fit("b", 0.1, mode = "tied"))
  expect_error(hierLangevin(fl, silent = TRUE), "must use the same lambda mode")
  fl <- list(a = make_lambda_fit("a", 0.05, mode = "fixed"), b = make_lambda_fit("b", 0.05))
  expect_error(hierLangevin(fl, silent = TRUE), "must use the same lambda mode")
  fl <- list(a = make_lambda_fit("a", 0.05, mode = "fixed"), b = make_lambda_fit("b", 0.1, mode = "fixed"))
  expect_error(hierLangevin(fl, silent = TRUE), "fitted with a fixed lambda must use the same value")
  # common fixed value and all tied are fine; lambda is not a Stage II parameter
  fl <- list(a = make_lambda_fit("a", 0.05, mode = "fixed"), b = make_lambda_fit("b", 0.05, mode = "fixed", beta = -1.3))
  h <- hierLangevin(fl, silent = TRUE)
  expect_equal(h$conditions$lambda, 0.05)
  expect_equal(h$conditions$lambda_mode, "fixed")
  expect_false(any(grepl("lambda", rownames(h$estimates$working))))
  fl <- list(a = make_lambda_fit("a", NA, mode = "tied"), b = make_lambda_fit("b", NA, mode = "tied", beta = -1.3))
  h <- hierLangevin(fl, silent = TRUE)
  expect_equal(h$conditions$lambda, Inf)
  expect_equal(h$conditions$lambda_mode, "tied")
})

test_that("hierLangevin stops if an individual's log_lambda SE is not finite", {
  fl <- list(a = make_lambda_fit("a", 0.05), b = make_lambda_fit("b", 0.1))
  fl$b$covariance$working["log_lambda", "log_lambda"] <- NaN
  expect_error(hierLangevin(fl, silent = TRUE), "standard error of log_lambda is not finite for individual 'b'")
})

test_that("getUD and regionProb integrate over the pooled lambda of a hierLangevin fit", {
  skip_on_cran()
  set.seed(1)
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30)
  cov1 <- terra::init(r, "y") / 30; names(cov1) <- "cov1"
  # signed distance to a vertical coast at x = 10 (restricted side: x < 10)
  bar <- terra::init(r, "x") - 10; names(bar) <- "coast"; attr(bar, "barLangevin") <- TRUE
  covs <- list(cov1 = cov1, coast = bar)
  lam <- c(0.005, 0.02, 0.08)
  fl <- setNames(lapply(1:3, function(i) make_lambda_fit(paste0("ind_", i), lam[i], beta = c(-1, -1.2, -0.8)[i], lse = 0.05)), paste0("ind_", 1:3))
  fl <- lapply(fl, function(f) { f$signatures <- list(covs = langevinSSM:::get_covs_signature(covs)); f })
  h <- hierLangevin(fl, silent = TRUE)
  ud <- suppressMessages(getUD(covs, h, individual = "all", nSims = 20, log = FALSE, plot = FALSE, show_progress = FALSE))
  land <- terra::ifel(bar <= 0, 1, 0)
  use_land <- sapply(paste0("UD_ID_ind_", 1:3), function(nm) terra::global(ud[[nm]] * land, "sum", na.rm = TRUE)[[1]])
  # weaker barrier (smaller lambda) -> more use of the restricted side
  expect_true(all(diff(use_land) < 0))
  # individual UD equals the closed form with that individual's lambda
  b1 <- h$estimates$random$beta_cov1$est$est[1]; l1 <- h$estimates$random$lambda$est$est[1]
  w <- terra::values(cov1)[, 1] * b1 - 0.5 * l1 * pmin(terra::values(bar)[, 1], 0)^2
  p <- exp(w - max(w)); p <- p / sum(p)
  expect_equal(as.numeric(terra::values(ud[["UD_ID_ind_1"]])[, 1]), p, tolerance = 1e-8)
  # population UD averages over log(lambda) ~ N(mu, sd^2): use of land lies inside the individual range
  pop_land <- terra::global(ud[["UD_Population"]] * land, "sum", na.rm = TRUE)[[1]]
  expect_true(pop_land > min(use_land) && pop_land < max(use_land))
  rp <- suppressMessages(regionProb(h, covs, mask = land, nSims = 50, individual = "ind_1", show_progress = FALSE))
  expect_equal(rp[["ID_ind_1"]]$Point_Estimate, use_land[1], tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(is.finite(rp[["Population"]]$Point_Estimate))
  # manual lambda = "estimate" without a pooled fit is rejected
  expect_error(getUD(covs, beta = 1, barrier = "coast", lambda = "estimate", plot = FALSE), "only supported for hierLangevin")
})
