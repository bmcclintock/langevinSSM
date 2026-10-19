# tests/testthat/test_obs_groups.R
#
# Observation-model extensions, all of which reduce to the original model at their defaults:
#   * obs_dist by error type ('epar' = error ellipse, 'sderr' = x.err/y.err) or location class
#   * tau by location class (observations with x.err/y.err)
#   * df / k by error type or location class
#   * error ellipse power law: smaj* = psi_maj M0 (smaj/M0)^b_maj, smin* = psi m0 (smin/m0)^b_min

grp_covs <- function() {
  r <- terra::rast(nrows = 30, ncols = 30, ext = c(-30, 30, -30, 30))
  set.seed(11)
  terra::values(r) <- stats::rnorm(900)
  list(cov = r)
}

# true track (overdamped, no measurement error)
grp_track <- function(n = 40, seed = 1) {
  covs <- grp_covs()
  set.seed(seed)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = covs, par = list(beta = 0.1, sigma = 1),
                                      timeStep = 0.1, obsPerAnimal = n, initialPosition = matrix(0, 1, 2)))
  list(sim = sim, covs = covs)
}

# LS/GPS data with location classes: alternating G / Argos classes
grp_ls <- function(n = 40, seed = 1, lcs = c("G", "3", "2", "1", "0", "A", "B")) {
  tr <- grp_track(n, seed)
  sim <- tr$sim
  sim$lc <- factor(rep(lcs, length.out = nrow(sim)))
  sim$x.err <- ifelse(sim$lc == "G", 0.05, 0.4); sim$y.err <- ifelse(sim$lc == "G", 0.05, 0.3)
  sim$smaj <- NA; sim$smin <- NA; sim$eor <- NA
  sim <- suppressMessages(addMeasurementError(sim))
  list(data = sim, covs = tr$covs)
}

# mixed data: GPS (LS, lc = G) and Argos error ellipses (KF, lc = 3, 2, ...)
grp_mixed <- function(n = 40, seed = 2, lc = TRUE) {
  tr <- grp_track(n, seed)
  sim <- tr$sim
  kf <- seq_len(nrow(sim)) %% 2 == 0
  set.seed(seed)
  sim$smaj <- ifelse(kf, stats::runif(nrow(sim), 0.2, 2), NA)
  sim$smin <- ifelse(kf, sim$smaj * stats::runif(nrow(sim), 0.1, 0.9), NA)
  sim$eor <- ifelse(kf, stats::runif(nrow(sim), 0, pi), NA)
  sim$x.err <- ifelse(kf, NA, 0.05); sim$y.err <- ifelse(kf, NA, 0.05)
  if (lc) sim$lc <- factor(ifelse(kf, c("3", "2", "1", "A", "B")[seq_len(nrow(sim)) %% 5 + 1], "G"))
  sim <- suppressMessages(addMeasurementError(sim))
  list(data = sim, covs = tr$covs)
}

# KF only
grp_kf <- function(n = 40, seed = 3) {
  tr <- grp_track(n, seed)
  set.seed(seed)
  sim <- suppressMessages(addMeasurementError(tr$sim, measurementError = list(smaj.sd = 1, smin.sd = 0.3)))
  list(data = sim, covs = tr$covs)
}

gfit <- function(d, par = list(), map = list(), obs_dist = "Gaussian", calcSE = FALSE, ...) {
  p <- c(list(beta = 0.1, sigma = 1), par)
  m <- c(list(beta = factor(NA), sigma = factor(NA)), map)
  suppressWarnings(suppressMessages(fitLangevin(d$data, model = "overdamped", obs_dist = obs_dist, spatialCovs = d$covs,
                                                par = p, map = m, calcSE = calcSE, silent = TRUE, ...)))
}

test_that("grouped parameters at the shared values reproduce the default likelihood", {
  d <- grp_ls()
  f0 <- gfit(d, list(tau = c(1.3, 0.8)), list(tau = factor(c(NA, NA))))
  lc <- levels(d$data$lc)
  tm <- matrix(rep(c(1.3, 0.8), length(lc)), 2, dimnames = list(NULL, lc))
  f1 <- gfit(d, list(tau = tm), list(tau = factor(rep(NA, length(tm)))))
  expect_equal(f1$objective, f0$objective, tolerance = 1e-10)
  # named list form, with "sderr" for all other classes
  f2 <- gfit(d, list(tau = list(G = c(1.3, 0.8), sderr = c(1.3, 0.8))), list(tau = factor(rep(NA, 4))))
  expect_equal(f2$objective, f0$objective, tolerance = 1e-10)
  expect_equal(f2$conditions$obs_spec$tau_labels, c("G", "sderr"))

  # t and Huber: df / k by group at the shared value
  for (od in c("t", "Huber")) {
    nm <- if (od == "t") "df" else "k"
    p0 <- list(tau = c(1, 1)); p0[[nm]] <- 3
    m0 <- list(tau = factor(c(NA, NA))); m0[[nm]] <- factor(NA)
    g0 <- gfit(d, p0, m0, obs_dist = od)
    p1 <- p0; p1[[nm]] <- c(G = 3, sderr = 3)
    m1 <- m0; m1[[nm]] <- factor(c(NA, NA))
    g1 <- gfit(d, p1, m1, obs_dist = od)
    expect_equal(g1$objective, g0$objective, tolerance = 1e-10)
  }

  # error ellipse power law at psi_maj = b_maj = b_min = 1
  k <- grp_kf()
  h0 <- gfit(k, list(psi = 1.4), list(psi = factor(NA)))
  h1 <- gfit(k, list(psi = 1.4, psi_maj = 1, b_maj = 1, b_min = 1),
             list(psi = factor(NA), psi_maj = factor(NA), b_maj = factor(NA), b_min = factor(NA)))
  expect_true(h1$conditions$obs_spec$ell_power)
  expect_equal(h1$objective, h0$objective, tolerance = 1e-10)
})

test_that("tau by location class enters the likelihood as tau of each observation's class", {
  d <- grp_ls(lcs = c("G", "A"))
  tauG <- c(0.7, 1.6); tauA <- c(2.1, 0.9)
  f <- gfit(d, list(tau = list(G = tauG, A = tauA)), list(tau = factor(rep(NA, 4))))
  # equivalent model: scale x.err / y.err directly and use tau = 1
  d2 <- d
  g <- d2$data$lc == "G"
  d2$data$x.err <- d2$data$x.err * ifelse(g, tauG[1], tauA[1])
  d2$data$y.err <- d2$data$y.err * ifelse(g, tauG[2], tauA[2])
  f2 <- gfit(d2, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))))
  expect_equal(f$objective, f2$objective, tolerance = 1e-8)
  expect_equal(rownames(f$estimates$natural)[grepl("^tau", rownames(f$estimates$natural))],
               c("tau_G_1", "tau_G_2", "tau_A_1", "tau_A_2"))
})

test_that("tau by location class is estimated and returned by getPar", {
  d <- grp_ls(n = 60, lcs = c("G", "A"))
  f <- gfit(d, list(tau = list(G = 1, A = 1)), list(tau = factor(c(1, 1, 2, 2))))
  expect_equal(f$convergence, 0)
  expect_equal(sum(names(f$par) == "l_tau"), 2)
  gp <- getPar(f)
  expect_true(is.matrix(gp$tau))
  expect_equal(colnames(gp$tau), c("G", "A"))
  expect_equal(gp$tau["x", ], gp$tau["y", ])
  # refit from getPar is at the optimum
  f2 <- gfit(d, list(tau = gp$tau), list(tau = factor(c(1, 1, 2, 2))))
  expect_equal(f2$objective, f$objective, tolerance = 1e-6)
})

test_that("df / k by error type works without location classes", {
  d <- grp_mixed(lc = FALSE)
  expect_true(all(is.na(d$data$lc)))
  f <- gfit(d, list(tau = c(1, 1), df = c(sderr = 20, epar = 4)), list(tau = factor(c(NA, NA)), df = factor(c(NA, NA))), obs_dist = "t")
  # direct check: the observation part of the likelihood uses df of each observation's type
  obj <- .make_obj(f, d$data, d$covs)
  dat <- obj$env$data
  expect_equal(dat$df_group, ifelse(is.na(d$data$smaj), 0L, 1L))
  expect_equal(dat$dist_obs, rep(1L, nrow(d$data)))
  expect_equal(rownames(f$estimates$natural)[grepl("^df", rownames(f$estimates$natural))], c("df_sderr", "df_epar"))
  # location-class names need lc
  expect_error(gfit(d, list(df = c(G = 5, epar = 5)), obs_dist = "t"), "no location classes")
})

test_that("df by location class: t likelihood equals the sum of the per-group likelihoods", {
  d <- grp_mixed()
  fG <- gfit(d, list(tau = c(1, 1), df = c(G = 30, epar = 3)), list(tau = factor(c(NA, NA)), df = factor(c(NA, NA))), obs_dist = "t")
  # with a fixed latent path the observation terms are separable: compare against two shared-df fits
  objG <- .make_obj(fG, d$data, d$covs)
  full <- objG$env$last.par.best   # only the latent path is free (all parameters fixed)
  # TMB tapes the data, so observation subsets need new objects: f_keep = process + observations in 'keep'
  mk <- function(f, keep) {
    dat <- .fit_tmb_data(f, d$data, d$covs)
    expect_equal(length(dat$isd), length(keep))
    dat$isd[!keep] <- 0L
    o <- TMB::MakeADFun(c(model = "langevinSSM", dat), f$tmb_setup$parList, map = f$tmb_setup$map,
                        random = f$tmb_setup$random, DLL = "langevinSSM_TMBExports", silent = TRUE)
    as.numeric(o$env$f(full))
  }
  sh <- function(df) gfit(d, list(tau = c(1, 1), df = df), list(tau = factor(c(NA, NA)), df = factor(NA)), obs_dist = "t")
  f30 <- sh(30); f3 <- sh(3)
  isG <- d$data$lc == "G"
  none <- rep(FALSE, nrow(d$data))
  ref <- mk(f30, isG) + mk(f3, !isG) - mk(f30, none)
  expect_equal(as.numeric(objG$env$f(full)), as.numeric(ref), tolerance = 1e-8)
})

test_that("error ellipse power law scales the axes as specified", {
  d <- grp_kf()
  pm <- 1.3; bM <- 0.5; bm <- 0.7; ps <- 2
  f <- gfit(d, list(psi = ps, psi_maj = pm, b_maj = bM, b_min = bm),
            list(psi = factor(NA), psi_maj = factor(NA), b_maj = factor(NA), b_min = factor(NA)))
  ref <- f$conditions$obs_spec$ell_ref
  expect_equal(unname(ref), unname(c(exp(mean(log(d$data$smaj))), exp(mean(log(d$data$smin))))))
  d2 <- d
  d2$data$smaj <- pm * ref[["maj"]] * (d$data$smaj / ref[["maj"]])^bM
  d2$data$smin <- ps * ref[["min"]] * (d$data$smin / ref[["min"]])^bm
  f2 <- gfit(d2)
  expect_equal(f$objective, f2$objective, tolerance = 1e-8)
  # invariant to scaleFactor
  f3 <- gfit(d, list(psi = ps, psi_maj = pm, b_maj = bM, b_min = bm),
             list(psi = factor(NA), psi_maj = factor(NA), b_maj = factor(NA), b_min = factor(NA)), scaleFactor = 10)
  expect_equal(getPar(f3)$mu, getPar(f)$mu, tolerance = 1e-5)
})

test_that("error ellipse power law parameters are estimated and reported", {
  d <- grp_kf(n = 80)
  f <- gfit(d, list(psi = 1, psi_maj = 1, b_maj = 1, b_min = 1), calcSE = TRUE)
  expect_true(all(c("l_psi_maj", "l_b_maj", "l_b_min") %in% names(f$par)))
  expect_true(all(c("psi_maj", "b_maj", "b_min") %in% rownames(f$estimates$natural)))
  gp <- getPar(f)
  expect_true(all(c("psi_maj", "b_maj", "b_min") %in% names(gp)))
  # default fit reports none of them
  f0 <- gfit(d)
  expect_false(any(c("psi_maj", "b_maj", "b_min") %in% rownames(f0$estimates$natural)))
  expect_false(any(c("l_psi_maj", "l_b_maj", "l_b_min") %in% names(f0$par)))
})

test_that("invalid or inapplicable group specifications are rejected", {
  d <- grp_ls(lcs = c("G", "A"))
  # missing class
  expect_error(gfit(d, list(tau = list(G = 1))), "no entry for .*location class A")
  # bad name
  expect_error(gfit(d, list(tau = list(G = 1, Q = 1))), "Invalid name")
  # tau only applies to x.err/y.err observations, so 'epar' is not a valid name
  expect_error(gfit(d, list(tau = list(epar = 1))), "Invalid name")
  # old development labels are no longer accepted
  expect_error(gfit(d, list(df = c(LS = 5, G = 5)), obs_dist = "t"), "Invalid name")
  expect_error(gfit(d, list(k = c(Argos = 5, G = 5)), obs_dist = "Huber"), "Invalid name")
  # map length
  expect_error(gfit(d, list(tau = list(G = 1, A = 1)), list(tau = factor(1:2))), "map\\$tau should be of length 4")
  # estimated class without observations warns
  expect_warning(suppressMessages(fitLangevin(d$data, model = "overdamped", spatialCovs = d$covs,
                                              par = list(beta = 0.1, sigma = 1, tau = list(G = 1, A = 1, B = 1)),
                                              map = list(beta = factor(NA), sigma = factor(NA)), calcSE = FALSE)),
                 "not identifiable")
  # power law without ellipses
  expect_error(gfit(d, list(b_maj = 1)), "error ellipse")
  # tau without LS data / psi_maj without KF data
  k <- grp_kf()
  expect_error(gfit(k, list(tau = list(G = 1))), "tau")
  # no lc column at all
  nolc <- grp_mixed(lc = FALSE)
  expect_error(gfit(nolc, list(tau = list(G = 1))), "no location classes")
})

test_that("grouped fits work with getPar, predLangevin, simLangevin and residuals", {
  d <- grp_mixed(n = 40)
  f <- gfit(d, list(tau = list(G = 1), k = c(G = 3, epar = 2), psi = 1, b_maj = 1),
            list(tau = factor(c(NA, NA))), obs_dist = "Huber", calcSE = TRUE)
  expect_equal(f$convergence, 0)
  gp <- getPar(f)
  expect_equal(names(gp$k), c("G", "epar"))
  expect_true(all(c("k_G", "k_epar", "tau_G_1", "tau_G_2") %in% rownames(f$estimates$natural)))
  expect_true("b_maj" %in% names(gp))
  # predLangevin keeps the groups and the fixed parameters
  pr <- suppressWarnings(suppressMessages(predLangevin(f, data = d$data, spatialCovs = d$covs, calcSE = FALSE)))
  expect_equal(as.numeric(pr$objective), as.numeric(f$objective), tolerance = 1e-6)
  expect_equal(pr$conditions$obs_spec$k_labels, c("G", "epar"))
  # simLangevin (conditional) uses the per-group values
  set.seed(1)
  s <- suppressWarnings(suppressMessages(simLangevin(f, data = d$data, spatialCovs = d$covs, conditional = TRUE)))
  expect_equal(nrow(s), nrow(d$data))
  expect_true(all(is.finite(s$x)))
  # residuals
  skip_on_cran()
  res <- suppressWarnings(suppressMessages(residuals(f, data = d$data, spatialCovs = d$covs, run_tests = FALSE)))
  expect_s3_class(res, "resLangevin")
  expect_gt(mean(is.finite(res$residual.x)), 0.8)
})

test_that("addMeasurementError and simLangevin accept grouped parameters", {
  tr <- grp_track(n = 3000)
  emf <- data.frame(lc = c("G", "A"), emf.x = c(0.1, 1), emf.y = c(0.1, 1), prob = c(0.5, 0.5))
  set.seed(3)
  s <- suppressMessages(addMeasurementError(tr$sim, par = list(tau = list(G = 1, A = 3)), measurementError = emf))
  rx <- (s$x - s$mu.x) / s$x.err
  expect_equal(stats::sd(rx[s$lc == "G"]), 1, tolerance = 0.1)
  expect_equal(stats::sd(rx[s$lc == "A"]), 3, tolerance = 0.3)
  # default tau reproduces the original simulation exactly
  set.seed(4); s1 <- suppressMessages(addMeasurementError(tr$sim, par = list(tau = c(1, 1)), measurementError = emf))
  set.seed(4); s2 <- suppressMessages(addMeasurementError(tr$sim, par = list(tau = list(G = 1, A = 1)), measurementError = emf))
  expect_equal(s1$x, s2$x)
  # power law at 1 = original ellipse simulation
  set.seed(5); k1 <- suppressMessages(addMeasurementError(tr$sim, par = list(psi = 1.5), measurementError = list(smaj.sd = 1, smin.sd = 0.3)))
  set.seed(5); k2 <- suppressMessages(addMeasurementError(tr$sim, par = list(psi = 1.5, b_maj = 1, b_min = 1, psi_maj = 1), measurementError = list(smaj.sd = 1, smin.sd = 0.3)))
  expect_equal(k1$x, k2$x)
  # df by error type with error ellipse simulation
  set.seed(6)
  k3 <- suppressMessages(addMeasurementError(tr$sim, par = list(df = c(epar = 4)), measurementError = list(smaj.sd = 1, smin.sd = 0.3), obs_dist = "t"))
  expect_true(all(is.finite(k3$x)))
  expect_error(suppressMessages(addMeasurementError(tr$sim, par = list(df = c(sderr = 4)), measurementError = list(smaj.sd = 1, smin.sd = 0.3), obs_dist = "t")),
               "no matching entry")
  # simLangevin with tau by class
  set.seed(7)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = tr$covs, par = list(beta = 0.1, sigma = 1, tau = list(G = 1, A = 2)),
                                      timeStep = 0.1, obsPerAnimal = 30, initialPosition = matrix(0, 1, 2), measurementError = emf))
  expect_true(all(c("G", "A") %in% sim$lc))
})

test_that("fits with tail_labels or Argos groups from development version r29 are upgraded", {
  spec <- .cond_obs_spec(list(obs_dist = "Huber", obs_spec = list(tau_labels = NULL, tail_labels = c("G", "Argos"))))
  expect_equal(spec$k_labels, c("G", "Argos"))
  expect_null(spec$df_labels)
  spec <- .cond_obs_spec(list(obs_dist = "t", obs_spec = list(tail_labels = c("LS", "KF"))))
  expect_equal(spec$df_labels, c("sderr", "epar"))
  expect_equal(.group_index(c("G", "Argos"), c("G", "A", "3", NA), c("sderr", "epar", "sderr", "sderr")), c(1L, 2L, 2L, NA))
})

# ---- obs_dist by error type / location class ---------------------------------------------------

test_that("obs_dist with the same distribution for every group reproduces the single obs_dist", {
  d <- grp_mixed()
  for (od in c("Gaussian", "t", "Huber")) {
    f0 <- gfit(d, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))), obs_dist = od)
    f1 <- gfit(d, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))), obs_dist = c(epar = od, sderr = od))
    f2 <- gfit(d, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))), obs_dist = list(G = od, epar = od))
    expect_equal(f1$objective, f0$objective, tolerance = 1e-8)
    expect_equal(f2$objective, f0$objective, tolerance = 1e-8)
  }
  # the formal default (all three choices) is still Gaussian
  expect_equal(.norm_obs_dist(c("Gaussian", "t", "Huber")), "Gaussian")
  expect_equal(.obs_dist_flag("Gaussian"), 0L)
})

test_that("mixed obs_dist: likelihood is the sum of the Gaussian and t parts", {
  d <- grp_mixed()
  fm <- gfit(d, list(tau = c(1, 1), df = 3), list(tau = factor(c(NA, NA)), df = factor(NA)),
             obs_dist = c(sderr = "Gaussian", epar = "t"))
  objm <- .make_obj(fm, d$data, d$covs)
  expect_equal(objm$env$data$dist_obs, ifelse(is.na(d$data$smaj), 0L, 1L))
  expect_equal(objm$env$data$obs_dist, 1L)
  full <- objm$env$last.par.best
  mk <- function(f, keep) {
    dat <- .fit_tmb_data(f, d$data, d$covs)
    dat$isd[!keep] <- 0L
    o <- TMB::MakeADFun(c(model = "langevinSSM", dat), f$tmb_setup$parList, map = f$tmb_setup$map,
                        random = f$tmb_setup$random, DLL = "langevinSSM_TMBExports", silent = TRUE)
    as.numeric(o$env$f(full))
  }
  fg <- gfit(d, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))), obs_dist = "Gaussian")
  ft <- gfit(d, list(tau = c(1, 1), df = 3), list(tau = factor(c(NA, NA)), df = factor(NA)), obs_dist = "t")
  ep <- !is.na(d$data$smaj)
  none <- rep(FALSE, nrow(d$data))
  ref <- mk(fg, !ep) + mk(ft, ep) - mk(fg, none)
  expect_equal(as.numeric(objm$env$f(full)), ref, tolerance = 1e-8)
})

test_that("mixed obs_dist estimates and reports only the parameters it uses", {
  d <- grp_mixed(n = 60)
  f <- gfit(d, list(tau = c(1, 1)), list(tau = factor(c(NA, NA))), obs_dist = c(epar = "Huber", sderr = "Gaussian"), calcSE = TRUE)
  expect_equal(f$convergence, 0)
  rn <- rownames(f$estimates$natural)
  expect_true("k" %in% rn)
  expect_false(any(grepl("^df", rn)))
  expect_equal(f$conditions$obs_dist, c(epar = "Huber", sderr = "Gaussian"))
  expect_output(print(f), "epar = Huber")
  # df cannot be given when no observations use t
  expect_error(gfit(d, list(df = 4), obs_dist = c(epar = "Huber", sderr = "Gaussian")), "unless obs_dist = 't'")
})

test_that("a location class entry in obs_dist overrides the error type entry", {
  d <- grp_mixed()
  f <- gfit(d, list(tau = c(1, 1), df = 3), list(tau = factor(c(NA, NA)), df = factor(NA)),
            obs_dist = c(sderr = "Gaussian", epar = "Gaussian", B = "t"))
  dat <- .make_obj(f, d$data, d$covs)$env$data
  expect_equal(dat$dist_obs, ifelse(d$data$lc == "B", 1L, 0L))
  # df and k grouped within the observations that use them
  g <- gfit(d, list(tau = c(1, 1), df = c(B = 3, epar = 10), k = 2),
            list(tau = factor(c(NA, NA)), df = factor(c(NA, NA)), k = factor(NA)),
            obs_dist = c(G = "Huber", epar = "t"))
  dat <- .make_obj(g, d$data, d$covs)$env$data
  isG <- d$data$lc == "G"
  expect_equal(dat$dist_obs, ifelse(isG, 2L, 1L))
  expect_equal(dat$obs_dist, 3L)
  expect_equal(dat$df_group[!isG], ifelse(d$data$lc[!isG] == "B", 0L, 1L))
  expect_equal(g$conditions$obs_spec$df_labels, c("B", "epar"))
  expect_equal(g$conditions$obs_spec$k_labels, NULL)
})

test_that("invalid obs_dist specifications are rejected", {
  d <- grp_mixed()
  expect_error(gfit(d, obs_dist = c(epar = "Huber")), "no entry for")
  expect_error(gfit(d, obs_dist = c(epar = "Cauchy", sderr = "t")), "Invalid distribution")
  expect_error(gfit(d, obs_dist = c(GPS = "t", epar = "t")), "Invalid name")
  expect_error(gfit(d, obs_dist = c("t", "Huber")), "obs_dist")
  nolc <- grp_mixed(lc = FALSE)
  expect_error(gfit(nolc, obs_dist = c(G = "t", epar = "Huber")), "no location classes")
})

test_that("addMeasurementError and simLangevin accept obs_dist by group", {
  tr <- grp_track(n = 3000)
  emf <- data.frame(lc = c("G", "A"), emf.x = c(0.1, 1), emf.y = c(0.1, 1), prob = c(0.5, 0.5))
  # same distribution for all groups = single obs_dist
  set.seed(8); s1 <- suppressMessages(addMeasurementError(tr$sim, par = list(df = 4), measurementError = emf, obs_dist = "t"))
  set.seed(8); s2 <- suppressMessages(addMeasurementError(tr$sim, par = list(df = 4), measurementError = emf, obs_dist = c(sderr = "t")))
  expect_equal(s1$x, s2$x)
  # Gaussian for G, t3 for A
  set.seed(9)
  s <- suppressMessages(addMeasurementError(tr$sim, par = list(df = 3), measurementError = emf, obs_dist = c(G = "Gaussian", A = "t")))
  rx <- (s$x - s$mu.x) / s$x.err
  r <- function(z) stats::sd(z) / stats::mad(z)
  expect_lt(abs(r(rx[s$lc == "G"]) - 1), 0.1)
  expect_gt(r(rx[s$lc == "A"]), 1.2)
  # simLangevin
  set.seed(10)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = tr$covs, par = list(beta = 0.1, sigma = 1, k = 1.5),
                                      timeStep = 0.1, obsPerAnimal = 30, initialPosition = matrix(0, 1, 2), measurementError = emf,
                                      obs_dist = c(G = "Gaussian", A = "Huber")))
  expect_true(all(is.finite(sim$x)))
})

test_that("mixed obs_dist fits work with predLangevin, simLangevin and residuals", {
  d <- grp_mixed(n = 40)
  od <- c(sderr = "Gaussian", epar = "Huber")
  f <- gfit(d, list(tau = c(1, 1), k = c("3" = 3, epar = 2)), list(tau = factor(c(NA, NA))), obs_dist = od)
  expect_equal(f$conditions$obs_spec$k_labels, c("3", "epar"))
  pr <- suppressWarnings(suppressMessages(predLangevin(f, data = d$data, spatialCovs = d$covs, calcSE = FALSE)))
  expect_equal(as.numeric(pr$objective), as.numeric(f$objective), tolerance = 1e-6)
  expect_equal(pr$conditions$obs_dist, od)
  set.seed(1)
  s <- suppressWarnings(suppressMessages(simLangevin(f, data = d$data, spatialCovs = d$covs, conditional = TRUE)))
  expect_true(all(is.finite(s$x)))
  skip_on_cran()
  res <- suppressWarnings(suppressMessages(residuals(f, data = d$data, spatialCovs = d$covs, run_tests = FALSE)))
  expect_s3_class(res, "resLangevin")
})
