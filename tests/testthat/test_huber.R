# tests/testthat/test_huber.R
#
# Pseudo-Huber (bivariate hyperbolic) observation model, obs_dist = "Huber":
#   log f(r) = -(k sqrt(k^2 + d2) - k^2) - log(2 pi) - 1/2 log|Sigma| - log(1 + 1/k^2),   d2 = r' Sigma^-1 r

ldhub <- function(r, S, k) {
  d2 <- drop(t(r) %*% solve(S) %*% r)
  -(k * sqrt(k^2 + d2) - k^2) - log(2 * pi) - 0.5 * log(det(S)) - log(1 + 1 / k^2)
}
ldgau <- function(r, S) {
  d2 <- drop(t(r) %*% solve(S) %*% r)
  -0.5 * d2 - log(2 * pi) - 0.5 * log(det(S))
}

huber_sim <- function(n = 30, k = 1.3, rho = 0.3, seed = 1) {
  set.seed(seed)
  r <- terra::rast(nrows = 30, ncols = 30, ext = c(-30, 30, -30, 30))
  terra::values(r) <- stats::rnorm(900)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = list(cov = r),
                                      par = list(beta = 0.1, sigma = 1), timeStep = 0.1,
                                      obsPerAnimal = n, initialPosition = matrix(0, 1, 2)))
  sim$x.err <- 0.5; sim$y.err <- 0.3
  sim <- suppressMessages(addMeasurementError(sim, par = list(tau = c(1, 1), rho_o = rho, k = k),
                                              obs_dist = "Huber"))
  list(data = sim, covs = list(cov = r))
}

# fit with everything but the latent states fixed and return the TMB object
fixed_obj <- function(hs, obs_dist, extra = list(), rho = 0.3) {
  p <- c(list(beta = 0.1, sigma = 1, tau = c(1, 1), rho_o = rho), extra)
  m <- list(beta = factor(NA), sigma = factor(NA), tau = factor(c(NA, NA)), rho_o = factor(NA))
  for (n in names(extra)) m[[n]] <- factor(NA)
  fit <- suppressWarnings(suppressMessages(
    fitLangevin(hs$data, model = "overdamped", obs_dist = obs_dist, spatialCovs = hs$covs, par = p, map = m,
                calcSE = FALSE, silent = TRUE)))
  list(fit = fit, obj = .make_obj(fit, hs$data, hs$covs))
}

test_that("the bivariate hyperbolic density is normalised", {
  k <- 1.7; S <- matrix(c(2, 0.6, 0.6, 1), 2)
  inner <- function(xx) stats::integrate(function(y) sapply(y, function(yy) exp(ldhub(c(xx, yy), S, k))), -Inf, Inf)$value
  I <- stats::integrate(function(x) sapply(x, inner), -Inf, Inf)$value
  expect_equal(I, 1, tolerance = 1e-5)
})

test_that("obs_dist = 'Huber' evaluates the pseudo-Huber likelihood exactly", {
  hs <- huber_sim()
  k <- 1.3; rho <- 0.3
  H <- fixed_obj(hs, "Huber", list(k = k))
  G <- fixed_obj(hs, "Gaussian")
  expect_equal(H$fit$conditions$obs_dist, "Huber")
  full <- H$obj$env$last.par.best
  mu <- t(matrix(full[H$obj$env$random], 2))
  S <- matrix(c(0.25, rho * 0.15, rho * 0.15, 0.09), 2)
  obs <- cbind(hs$data$x, hs$data$y)
  ok <- which(!is.na(obs[, 1]))
  ref <- -sum(sapply(ok, function(i) ldhub(obs[i, ] - mu[i, ], S, k))) +
          sum(sapply(ok, function(i) ldgau(obs[i, ] - mu[i, ], S)))
  # same process model and latent path, so the difference is the observation part only
  expect_equal(as.numeric(H$obj$env$f(full) - G$obj$env$f(full)), ref, tolerance = 1e-8)
  # k is reported on the natural scale
  expect_equal(H$fit$estimates$natural["k", "Estimate"], k)
})

test_that("obs_dist = 'Huber' tends to the Gaussian model as k grows", {
  hs <- huber_sim()
  G <- fixed_obj(hs, "Gaussian")
  H <- fixed_obj(hs, "Huber", list(k = 500))
  expect_equal(H$fit$objective, G$fit$objective, tolerance = 1e-3)
})

test_that("the pseudo-Huber observation term is convex in the latent locations (unlike the t)", {
  # with beta = 0 the process model is Gaussian, so the inner Hessian must be positive definite at any
  # latent path for Gaussian and Huber errors, while t errors can make it indefinite away from the mode
  set.seed(2)
  r <- terra::rast(nrows = 30, ncols = 30, ext = c(-30, 30, -30, 30)); terra::values(r) <- 0
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = list(cov = r),
                                      par = list(beta = 0, sigma = 1), timeStep = 0.1,
                                      obsPerAnimal = 30, initialPosition = matrix(0, 1, 2)))
  sim$x.err <- 0.5; sim$y.err <- 0.3
  sim <- suppressMessages(addMeasurementError(sim, par = list(tau = c(1, 1), rho_o = 0.3, k = 1),
                                              obs_dist = "Huber"))
  hs <- list(data = sim, covs = list(cov = r))
  mineig <- function(o) {
    set.seed(3); rr <- o$env$random; b <- o$env$last.par.best
    min(sapply(1:20, function(i) {
      f2 <- b; f2[rr] <- b[rr] + stats::rnorm(length(rr), 0, 5)
      min(eigen(as.matrix(o$env$spHess(f2, random = TRUE)), symmetric = TRUE, only.values = TRUE)$values)
    }))
  }
  fx <- function(od, extra) {
    p <- c(list(beta = 0, sigma = 1, tau = c(1, 1), rho_o = 0.3), extra)
    m <- list(beta = factor(NA), sigma = factor(NA), tau = factor(c(NA, NA)), rho_o = factor(NA))
    for (n in names(extra)) m[[n]] <- factor(NA)
    f <- suppressWarnings(suppressMessages(fitLangevin(sim, model = "overdamped", obs_dist = od, spatialCovs = hs$covs,
                                                       par = p, map = m, calcSE = FALSE, silent = TRUE)))
    .make_obj(f, sim, hs$covs)
  }
  expect_gt(mineig(fx("Huber", list(k = 1))), 0)
  expect_lt(mineig(fx("t", list(df = 3))), 0)
})

test_that("addMeasurementError(obs_dist = 'Huber') samples from the bivariate hyperbolic distribution", {
  set.seed(29)
  n <- 20000; k <- 1.2; rho <- 0.4
  d0 <- data.frame(id = 1, date = seq_len(n), dt = 1, mu.x = 0, mu.y = 0, x.err = 0.5, y.err = 0.3,
                   smaj = NA, smin = NA, eor = NA)
  class(d0) <- c("dataLangevin", "data.frame")
  e <- suppressMessages(addMeasurementError(d0, par = list(tau = c(1, 1), rho_o = rho, k = k), obs_dist = "Huber"))
  Z <- cbind(e$x, e$y)
  S <- matrix(c(0.25, rho * 0.15, rho * 0.15, 0.09), 2)
  p <- k^2 / (k^2 + 1)
  Es2 <- p * (1 + 2 / k^2 + 2 / k^4) + (1 - p) * (1 + 4 / k^2 + 6 / k^4)
  expect_equal(c(stats::cov(Z)), c(S * k^2 * (Es2 - 1) / 2), tolerance = 0.05)
  d <- sqrt(rowSums((Z %*% solve(S)) * Z))
  Fd <- function(d) { t <- sqrt(1 + d^2 / k^2) - 1; 1 - exp(-k^2 * t) * (p + (1 - p) * (1 + k^2 * t)) }
  expect_gt(suppressWarnings(stats::ks.test(d, Fd))$p.value, 0.001)
})

test_that("obs_dist = 'Huber' parameter handling", {
  hs <- huber_sim(n = 10)
  expect_error(suppressMessages(fitLangevin(hs$data, model = "overdamped", obs_dist = "Huber", spatialCovs = hs$covs,
                                            par = list(beta = 0, sigma = 1, k = -1))), "k")
  expect_error(suppressMessages(fitLangevin(hs$data, model = "overdamped", obs_dist = "Gaussian", spatialCovs = hs$covs,
                                            par = list(beta = 0, sigma = 1, k = 2))), "Huber")
  expect_error(suppressMessages(fitLangevin(hs$data, model = "overdamped", obs_dist = "t", spatialCovs = hs$covs,
                                            par = list(beta = 0, sigma = 1, k = 2))), "Huber")
  expect_error(suppressMessages(fitLangevin(hs$data, model = "overdamped", obs_dist = "Huber", spatialCovs = hs$covs,
                                            par = list(beta = 0, sigma = 1, df = 2))), "'t'")
  # k defaults to 2 and is fixed when it is not given
  f <- suppressWarnings(suppressMessages(fitLangevin(hs$data, model = "overdamped", obs_dist = "Huber",
                                                     spatialCovs = hs$covs, par = list(beta = 0, sigma = 1),
                                                     calcSE = FALSE, silent = TRUE)))
  expect_equal(f$estimates$natural["k", "Estimate"], 2)
  expect_true(is.na(f$tmb_setup$map$l_k))
  expect_null(getPar(f)$k)
})

test_that("obs_dist = 'Huber' estimates k and the observation parameters", {
  skip_on_cran()
  set.seed(5)
  r <- terra::rast(nrows = 40, ncols = 40, ext = c(-40, 40, -40, 40)); terra::values(r) <- 0
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = list(cov = r),
                                      par = list(beta = 0, sigma = 1), timeStep = 0.1,
                                      obsPerAnimal = 300, initialPosition = matrix(0, 1, 2)))
  sim$x.err <- 0.5; sim$y.err <- 0.3
  sim <- suppressMessages(addMeasurementError(sim, par = list(tau = c(1, 1), rho_o = 0.4, k = 1.2),
                                              obs_dist = "Huber"))
  f <- suppressWarnings(suppressMessages(fitLangevin(sim, model = "overdamped", obs_dist = "Huber",
                                                     spatialCovs = list(cov = r),
                                                     par = list(beta = 0, sigma = 1, tau = c(1, 1), rho_o = 0, k = 3),
                                                     map = list(beta = factor(NA)), silent = TRUE)))
  expect_equal(f$convergence, 0)
  est <- f$estimates$natural
  expect_true(all(is.finite(est[, "Std. Error"])))
  expect_lt(abs(est["k", "Estimate"] - 1.2), 3 * est["k", "Std. Error"])
  expect_lt(abs(est["rho_o", "Estimate"] - 0.4), 3 * est["rho_o", "Std. Error"])
  expect_true("k" %in% names(getPar(f)))
})

test_that("the pseudo-Huber one-step CDF matches numerical integration", {
  # x-marginal GH(3/2, k^2, k^2) and y|x-conditional GH(1, k^2 + zx^2, k^2), standardised
  ldgh1 <- function(z, lam, chi, psi) 0.5 * lam * log(psi / chi) - 0.5 * log(2 * pi) -
    (log(besselK(sqrt(chi * psi), lam, TRUE)) - sqrt(chi * psi)) + 0.5 * (lam - 0.5) * log((chi + z^2) / psi) +
    (log(besselK(sqrt(psi * (chi + z^2)), lam - 0.5, TRUE)) - sqrt(psi * (chi + z^2)))
  ref <- function(z, lam, chi, psi) stats::integrate(function(u) exp(ldgh1(u, lam, chi, psi)), z, Inf, rel.tol = 1e-12)$value
  hs <- huber_sim(n = 12)
  k <- 1.3; rho <- 0.3
  H <- fixed_obj(hs, "Huber", list(k = k))
  mu <- t(matrix(H$obj$env$last.par.best[H$obj$env$random], 2))
  S <- matrix(c(0.25, rho * 0.15, rho * 0.15, 0.09), 2)
  obs <- cbind(hs$data$x, hs$data$y); ok <- which(!is.na(obs[, 1]))
  zx <- (obs[ok, 1] - mu[ok, 1]) / sqrt(S[1, 1])
  vyx <- S[2, 2] - S[1, 2]^2 / S[1, 1]
  zy <- (obs[ok, 2] - mu[ok, 2] - S[1, 2] / S[1, 1] * (obs[ok, 1] - mu[ok, 1])) / sqrt(vyx)
  ux <- 1 - sapply(zx, function(z) ref(z, 1.5, k^2, k^2))
  uy <- 1 - mapply(function(z, a) ref(z, 1, k^2 + a^2, k^2), zy, zx)
  # with the latent path fixed (no random effects), oneStepPredict(method = "cdf") returns qnorm of these CDFs
  m <- H$fit$tmb_setup$map; m$mu <- factor(rep(NA, length(H$fit$tmb_setup$parList$mu)))
  dat <- .fit_tmb_data(H$fit, hs$data, hs$covs)
  pl <- H$fit$tmb_setup$parList; pl$mu <- matrix(H$obj$env$last.par.best[H$obj$env$random], 2)
  o <- TMB::MakeADFun(c(model = "langevinSSM", dat), pl, map = m, DLL = "langevinSSM_TMBExports", silent = TRUE)
  osa <- suppressWarnings(TMB::oneStepPredict(o, observation.name = "Y", data.term.indicator = "keep",
                                              method = "cdf", discrete = FALSE, trace = FALSE))
  r <- matrix(osa$residual, 2)[, ok]
  expect_equal(r[1, ], stats::qnorm(ux), tolerance = 1e-4)
  expect_equal(r[2, ], stats::qnorm(uy), tolerance = 1e-4)
  skip_on_cran()
  # residuals() accepts method = "cdf" for Huber fits
  res <- suppressWarnings(suppressMessages(residuals(H$fit, data = hs$data, spatialCovs = hs$covs, method = "cdf",
                                                     run_tests = FALSE)))
  expect_s3_class(res, "resLangevin")
  expect_gt(sum(is.finite(res$residual.x)), 0.8 * length(ok))
  expect_gt(sum(is.finite(res$residual.y)), 0.8 * length(ok))
})
