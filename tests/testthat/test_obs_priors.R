# tests/testthat/test_obs_priors.R
#
# Priors on the observation parameters: group names, tied elements counted once, and errors for fixed parameters.

pri_covs <- function() {
  r <- terra::rast(nrows = 30, ncols = 30, ext = c(-30, 30, -30, 30))
  set.seed(11)
  terra::values(r) <- stats::rnorm(900)
  list(cov = r)
}

# overdamped track with x.err/y.err observations of location classes G, 3 and B
pri_ls <- function(n = 40, seed = 1) {
  covs <- pri_covs()
  set.seed(seed)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = covs, par = list(beta = 0.1, sigma = 1),
                                      timeStep = 0.1, obsPerAnimal = n, initialPosition = matrix(0, 1, 2)))
  sim$lc <- factor(rep(c("G", "3", "B"), length.out = nrow(sim)))
  sim$x.err <- ifelse(sim$lc == "G", 0.05, 0.4); sim$y.err <- ifelse(sim$lc == "G", 0.05, 0.3)
  sim$smaj <- NA; sim$smin <- NA; sim$eor <- NA
  list(data = suppressMessages(addMeasurementError(sim)), covs = covs)
}

# GPS (x.err/y.err, class G) and Argos error ellipses (classes 3, 2, 1, A, B)
pri_mixed <- function(n = 40, seed = 2) {
  covs <- pri_covs()
  set.seed(seed)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = covs, par = list(beta = 0.1, sigma = 1),
                                      timeStep = 0.1, obsPerAnimal = n, initialPosition = matrix(0, 1, 2)))
  kf <- seq_len(nrow(sim)) %% 2 == 0
  sim$smaj <- ifelse(kf, stats::runif(nrow(sim), 0.2, 2), NA)
  sim$smin <- ifelse(kf, sim$smaj * stats::runif(nrow(sim), 0.1, 0.9), NA)
  sim$eor <- ifelse(kf, stats::runif(nrow(sim), 0, pi), NA)
  sim$x.err <- ifelse(kf, NA, 0.05); sim$y.err <- ifelse(kf, NA, 0.05)
  sim$lc <- factor(ifelse(kf, c("3", "2", "1", "A", "B")[seq_len(nrow(sim)) %% 5 + 1], "G"))
  list(data = suppressMessages(addMeasurementError(sim)), covs = covs)
}

# fit with beta and sigma fixed; iter.max = 0 evaluates the objective at the starting values
pfit <- function(d, par = list(), map = list(), prior = NULL, obs_dist = "Gaussian", eval_only = FALSE) {
  ctrl <- if (eval_only) list(iter.max = 0, eval.max = 1) else list(trace = 0, iter.max = 1000, eval.max = 1000)
  quiet <- function(w) if (grepl("converge|NaN|NA/", conditionMessage(w))) invokeRestart("muffleWarning")
  withCallingHandlers(suppressMessages(fitLangevin(d$data, model = "overdamped", obs_dist = obs_dist, spatialCovs = d$covs,
                               par = c(list(beta = 0.1, sigma = 1), par), map = c(list(beta = factor(NA), sigma = factor(NA)), map),
                               prior = prior, calcSE = FALSE, silent = TRUE, control = ctrl)), warning = quiet)
}
pr1 <- function(rn, m = 0.5, s = 0.1) data.frame(mean = rep(m, length(rn)), sd = rep(s, length(rn)), row.names = rn)
pri <- function(f, wn) f$tmb_setup$priors[[paste0("prior_mean_", wn)]]
tau3 <- list(B = 1, G = 1, "3" = 1)   # groups in this order: B, G, 3

test_that("tau priors can be given by group name and agree with the positional names", {
  d <- pri_ls()
  # columns of tau: B (elements 1, 2), G (3, 4), 3 (5, 6)
  expect_equal(pri(pfit(d, list(tau = tau3), prior = pr1("l_tau.G"), eval_only = TRUE), "l_tau"), c(NA, NA, 0.5, 0.5, NA, NA))
  expect_equal(pri(pfit(d, list(tau = tau3), prior = pr1("l_tau.x.3"), eval_only = TRUE), "l_tau"), c(NA, NA, NA, NA, 0.5, NA))
  expect_equal(pri(pfit(d, list(tau = tau3), prior = pr1("l_tau.y.B"), eval_only = TRUE), "l_tau"), c(NA, 0.5, NA, NA, NA, NA))
  expect_equal(pri(pfit(d, list(tau = tau3), prior = pr1("l_tau.x"), eval_only = TRUE), "l_tau"), c(0.5, NA, 0.5, NA, 0.5, NA))
  expect_equal(pri(pfit(d, list(tau = tau3), prior = pr1("l_tau_5"), eval_only = TRUE), "l_tau"), c(NA, NA, NA, NA, 0.5, NA))
  # same objective by name and by position
  a <- pfit(d, list(tau = tau3), prior = pr1("l_tau.x.3"), eval_only = TRUE)
  b <- pfit(d, list(tau = tau3), prior = pr1("l_tau_5"), eval_only = TRUE)
  expect_equal(a$objective, b$objective)
  # a later, more specific row overrides an earlier one
  p <- data.frame(mean = c(0, 1), sd = c(1, 0.1), row.names = c("l_tau", "l_tau.G"))
  expect_equal(pri(pfit(d, list(tau = tau3), prior = p, eval_only = TRUE), "l_tau"), c(0, 0, 1, 1, 0, 0))
  # shared tau (no groups): x and y names
  expect_equal(pri(pfit(d, list(tau = c(1, 1)), prior = pr1("l_tau.y"), eval_only = TRUE), "l_tau"), c(NA, 0.5))
  # unknown group
  expect_error(pfit(d, list(tau = tau3), prior = pr1("l_tau.A"), eval_only = TRUE), "Groups of observation parameters for this model")
})

test_that("a prior by group name pulls only that group", {
  d <- pri_ls()
  mp <- list(tau = factor(c(1, 1, 2, 2, 3, 3)))
  f0 <- pfit(d, list(tau = tau3), map = mp)
  f1 <- pfit(d, list(tau = tau3), map = mp, prior = pr1("l_tau.x.B", m = log(3), s = 0.01))
  t0 <- getPar(f0)$tau; t1 <- getPar(f1)$tau
  # B (x and y tied) is pinned at the prior mean; it was far from it without the prior
  expect_equal(unname(t1[, "B"]), c(3, 3), tolerance = 0.02)
  expect_true(abs(t0["x", "B"] - 3) > 1)
  expect_equal(pri(f1, "l_tau"), c(log(3), NA, NA, NA, NA, NA))
})

test_that("df and k priors can be given by group name", {
  d <- pri_mixed()
  f <- pfit(d, list(df = c(epar = 5, B = 3)), obs_dist = c(sderr = "Gaussian", epar = "t"), prior = pr1("l_df.B"), eval_only = TRUE)
  expect_equal(pri(f, "l_df"), c(NA, 0.5))
  f <- pfit(d, list(df = c(epar = 5, B = 3)), obs_dist = c(sderr = "Gaussian", epar = "t"), prior = pr1("l_df.epar"), eval_only = TRUE)
  expect_equal(pri(f, "l_df"), c(0.5, NA))
  f <- pfit(d, list(k = c(epar = 2, B = 1)), obs_dist = c(sderr = "Gaussian", epar = "Huber"), prior = pr1("l_k.B"), eval_only = TRUE)
  expect_equal(pri(f, "l_k"), c(NA, 0.5))
  # one tail parameter per group: no x/y names
  expect_error(pfit(d, list(df = c(epar = 5, B = 3)), obs_dist = c(sderr = "Gaussian", epar = "t"), prior = pr1("l_df.x"), eval_only = TRUE),
               "is invalid")
})

test_that("tied elements receive their prior once", {
  d <- pri_ls()
  tie <- list(tau = factor(c(1, 1, 2, 2, 3, 3)))
  g0 <- pfit(d, list(tau = tau3), map = tie, eval_only = TRUE)
  one <- -stats::dnorm(0, 1, 1, log = TRUE)
  # base name over 3 groups tied x = y: 3 estimated parameters, 3 prior terms (not 6)
  g1 <- pfit(d, list(tau = tau3), map = tie, prior = pr1("l_tau", m = 1, s = 1), eval_only = TRUE)
  expect_equal(g1$objective - g0$objective, 3 * one, tolerance = 1e-8)
  # group name covering a tied pair: one term
  g2 <- pfit(d, list(tau = tau3), map = tie, prior = pr1("l_tau.G", m = 1, s = 1), eval_only = TRUE)
  expect_equal(g2$objective - g0$objective, one, tolerance = 1e-8)
  # the same prior on both elements from different rows: warning, applied once
  expect_warning(g3 <- pfit(d, list(tau = tau3), map = tie, prior = pr1(c("l_tau_1", "l_tau_2"), m = 1, s = 1), eval_only = TRUE),
                 "tied by 'map'")
  expect_equal(g3$objective - g0$objective, one, tolerance = 1e-8)
  # ties across groups (one tau for B and 3)
  tie2 <- list(tau = factor(c(1, 2, 3, 4, 1, 2)))
  h0 <- pfit(d, list(tau = tau3), map = tie2, eval_only = TRUE)
  h1 <- pfit(d, list(tau = tau3), map = tie2, prior = pr1("l_tau", m = 1, s = 1), eval_only = TRUE)
  expect_equal(h1$objective - h0$objective, 4 * one, tolerance = 1e-8)
  # without ties every element is a parameter
  u0 <- pfit(d, list(tau = tau3), eval_only = TRUE)
  u1 <- pfit(d, list(tau = tau3), prior = pr1("l_tau", m = 1, s = 1), eval_only = TRUE)
  expect_equal(u1$objective - u0$objective, 6 * one, tolerance = 1e-8)
})

test_that("priors on fixed or unused parameters are rejected", {
  d <- pri_ls()
  # fixed by map
  expect_error(pfit(d, list(tau = tau3), map = list(tau = factor(rep(NA, 6))), prior = pr1("l_tau"), eval_only = TRUE),
               "not estimated")
  expect_error(pfit(d, list(tau = tau3), map = list(tau = factor(c(1, 2, NA, NA, 3, 4))), prior = pr1("l_tau.G"), eval_only = TRUE),
               "'l_tau.G' \\(of 'l_tau'\\), which is not estimated")
  expect_error(pfit(d, list(tau = tau3), map = list(tau = factor(c(1, 2, NA, NA, 3, 4))), prior = pr1("l_tau_3"), eval_only = TRUE),
               "not estimated")
  # base name with some fixed elements: only the estimated ones get the prior
  f <- pfit(d, list(tau = tau3), map = list(tau = factor(c(1, 2, NA, NA, 3, 4))), prior = pr1("l_tau"), eval_only = TRUE)
  expect_equal(pri(f, "l_tau"), c(0.5, 0.5, NA, NA, 0.5, 0.5))
  # not in par (fixed at the default)
  expect_error(pfit(d, prior = pr1("l_tau"), eval_only = TRUE), "'tau' is only estimated if it is included in 'par'")
  expect_error(pfit(d, prior = pr1("l_rho_o"), eval_only = TRUE), "not estimated")
  # not used by the model
  expect_error(pfit(d, prior = pr1("l_df"), eval_only = TRUE), "obs_dist = \"t\"")
  expect_error(pfit(d, prior = pr1("l_k"), eval_only = TRUE), "obs_dist = \"Huber\"")
  dm <- pri_mixed()
  expect_error(pfit(dm, list(psi = 1), prior = pr1("l_b_maj"), eval_only = TRUE), "power law")
  # movement parameters fixed by map
  expect_error(pfit(d, prior = pr1("log_sigma"), eval_only = TRUE), "'log_sigma', which is not estimated")
})

test_that("priors on the power law and psi work when estimated", {
  d <- pri_mixed()
  pl <- list(psi = 1, psi_maj = 1, b_maj = 1, b_min = 1)
  f0 <- pfit(d, pl, eval_only = TRUE)
  f1 <- pfit(d, pl, prior = pr1(c("l_psi_maj", "l_b_min"), m = 1, s = 1), eval_only = TRUE)
  expect_equal(f1$objective - f0$objective, -2 * stats::dnorm(0, 1, 1, log = TRUE), tolerance = 1e-8)
  # fixing one power law parameter by map makes its prior an error
  expect_error(pfit(d, pl, map = list(b_maj = factor(NA)), prior = pr1("l_b_maj"), eval_only = TRUE), "not estimated")
})

test_that("predLangevin and residuals work with observation priors", {
  d <- pri_ls()
  f <- pfit(d, list(tau = tau3), map = list(tau = factor(c(1, 1, 2, 2, 3, 3))), prior = pr1("l_tau.G", m = 0, s = 0.5))
  p <- suppressWarnings(suppressMessages(predLangevin(f, d$data, d$covs, calcSE = FALSE)))
  expect_s3_class(p, "predLangevin")
  r <- suppressWarnings(suppressMessages(residuals(f, d$data, d$covs, run_tests = FALSE)))
  expect_s3_class(r, "resLangevin")
})
