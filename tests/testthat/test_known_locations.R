# Observations with coordinates but no measurement error information are known locations (mu fixed at the
# observation), whether all or only some observations lack measurement error.

known_sim <- function(n = 60, seed = 1) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(-500, 500, -500, 500))
  terra::values(r) <- runif(terra::ncell(r))
  names(r) <- "habitat"
  covs <- list(habitat = r)
  dat <- suppressMessages(simLangevin(obsPerAnimal = n, model = "underdamped", spatialCovs = covs,
                                      par = list(beta = 0, sigma = 5, gamma = 0.5)))
  dat$date <- as.POSIXct("2024-01-01", tz = "UTC") + seq_len(nrow(dat)) * 3600
  dat$id <- "A"; dat$lc <- "G"
  class(dat) <- "data.frame"
  list(dat = dat, covs = covs)
}

fit_quiet <- function(...) suppressMessages(suppressWarnings(fitLangevin(...)))

test_that("mixed data: locations without measurement error are fixed at the observations", {
  s <- known_sim()
  d <- s$dat; d$x.err <- d$y.err <- 2
  k <- c(5, 17, 30, 44)
  d$x.err[k] <- d$y.err[k] <- NA
  fd <- formatData(d, time.unit = "hours")
  expect_warning(fit <- suppressMessages(fitLangevin(fd, model = "underdamped", spatialCovs = s$covs,
                                                     par = list(beta = 0, sigma = 5, gamma = 0.5))),
                 "4 of 60 observed locations have no measurement error information")
  kk <- which(is.na(fd$x.err) & !is.na(fd$x))
  expect_length(kk, 4)
  mu <- fit$estimates$random$mu
  expect_equal(mu$est$mu.x[kk], fd$x[kk])
  expect_equal(mu$est$mu.y[kk], fd$y[kk])
  expect_true(all(mu$se$mu.x[kk] == 0) && all(mu$se$mu.y[kk] == 0))
  # the other locations are still estimated
  oo <- setdiff(seq_len(nrow(fd)), kk)
  expect_true(all(mu$se$mu.x[oo] > 0))
  # the map fixes exactly the known locations; velocities stay random
  m <- matrix(as.character(fit$tmb_setup$map$mu), nrow = 2)
  expect_equal(which(is.na(m[1, ])), kk)
  expect_true(all(c("mu", "vel") %in% fit$tmb_setup$random))
})

test_that("mixed data: fixing a known location equals giving it a negligible error", {
  s <- known_sim(n = 40, seed = 3)
  d <- s$dat; d$x.err <- d$y.err <- 3
  k <- seq(2, 40, by = 4)
  dk <- d; dk$x.err[k] <- dk$y.err[k] <- NA
  dt <- d; dt$x.err[k] <- dt$y.err[k] <- 1e-5
  p <- list(beta = 0, sigma = 5, gamma = 0.5)
  fk <- fit_quiet(formatData(dk, time.unit = "hours"), spatialCovs = s$covs, par = p)
  ft <- fit_quiet(formatData(dt, time.unit = "hours"), spatialCovs = s$covs, par = p)
  expect_equal(fk$estimates$natural$Estimate, ft$estimates$natural$Estimate, tolerance = 1e-3)
  # With a negligible error the observation term of each such location (log(2*pi) + 2*log(1e-5)) cancels its
  # contribution to the Laplace approximation (0.5*log det H - log(2*pi)), so the marginal likelihoods coincide
  expect_equal(fk$objective, ft$objective, tolerance = 1e-4)
})

test_that("all observations without measurement error: unchanged (mu fixed, vel random)", {
  s <- known_sim(n = 50, seed = 2)
  d <- s$dat
  fd <- formatData(d, predTimes = data.frame(id = "A", date = d$date[10] + 1800), time.unit = "hours")
  expect_silent(fit <- suppressMessages(fitLangevin(fd, spatialCovs = s$covs, par = list(beta = 0, sigma = 5, gamma = 0.5))))
  obs <- !is.na(fd$x)
  expect_equal(fit$estimates$random$mu$est$mu.x[obs], fd$x[obs])
  m <- matrix(as.character(fit$tmb_setup$map$mu), nrow = 2)
  expect_equal(which(!is.na(m[1, ])), which(!obs))     # only the prediction time is estimated
  expect_true(all(c("mu", "vel") %in% fit$tmb_setup$random))
  # without prediction times, mu is not a random effect at all
  fit2 <- fit_quiet(formatData(d, time.unit = "hours"), spatialCovs = s$covs, par = list(beta = 0, sigma = 5, gamma = 0.5))
  expect_identical(fit2$tmb_setup$random, "vel")
})

test_that("overdamped mixed data: known locations fixed", {
  s <- known_sim(n = 50, seed = 4)
  d <- s$dat; d$x.err <- d$y.err <- 2; k <- c(3, 20, 35)
  d$x.err[k] <- d$y.err[k] <- NA
  fd <- formatData(d, time.unit = "hours")
  fit <- fit_quiet(fd, model = "overdamped", spatialCovs = s$covs, par = list(beta = 0, sigma = 5))
  kk <- which(is.na(fd$x.err))
  expect_equal(fit$estimates$random$mu$est$mu.x[kk], fd$x[kk])
  expect_identical(fit$tmb_setup$random, "mu")
})

test_that("duplicated times: an observation with error shares the fixed state of a known location", {
  s <- known_sim(n = 40, seed = 3)
  d <- s$dat; d$x.err <- d$y.err <- 2
  d$x.err[12] <- d$y.err[12] <- NA
  dup <- d[12, ]; dup$x <- dup$x + 1; dup$x.err <- dup$y.err <- 2     # same time, with error
  d <- rbind(d, dup); d <- d[order(d$date), ]
  fd <- suppressWarnings(formatData(d, time.unit = "hours"))
  expect_warning(fit <- suppressMessages(fitLangevin(fd, spatialCovs = s$covs, par = list(beta = 0, sigma = 5, gamma = 0.5))),
                 "1 observation\\(s\\) with measurement error at the same time")
  kk <- which(is.na(fd$x.err) & !is.na(fd$x)); j <- which(fd$date == fd$date[kk] & !is.na(fd$x.err))
  mu <- fit$estimates$random$mu$est
  expect_equal(mu$mu.x[c(kk, j)], rep(fd$x[kk], 2))
  expect_equal(mu$mu.y[c(kk, j)], rep(fd$y[kk], 2))
})

test_that("duplicated times: identical known locations are allowed, different ones are an error", {
  s <- known_sim(n = 30, seed = 6)
  d <- s$dat; d$x.err <- d$y.err <- 2
  d$x.err[8] <- d$y.err[8] <- NA
  d2 <- rbind(d, d[8, ]); d2 <- d2[order(d2$date), ]
  p <- list(beta = 0, sigma = 5, gamma = 0.5)
  fit <- fit_quiet(suppressWarnings(formatData(d2, time.unit = "hours")), spatialCovs = s$covs, par = p)
  expect_s3_class(fit, "fitLangevin")
  dup <- d[8, ]; dup$x <- dup$x + 1
  d3 <- rbind(d, dup); d3 <- d3[order(d3$date), ]
  expect_error(suppressWarnings(fitLangevin(suppressWarnings(formatData(d3, time.unit = "hours")), spatialCovs = s$covs, par = p)),
               "have no measurement error information, so they are known locations, but their coordinates differ")
})

test_that("priors cannot be put on known locations, but can on estimated ones", {
  s <- known_sim(n = 30, seed = 7)
  d <- s$dat; d$x.err <- d$y.err <- 2; d$x.err[10] <- d$y.err[10] <- NA
  fd <- formatData(d, time.unit = "hours")
  p <- list(beta = 0, sigma = 5, gamma = 0.5)
  expect_error(suppressWarnings(fitLangevin(fd, spatialCovs = s$covs, par = p,
                                            prior = data.frame(mean = fd$x[10], sd = 1, row.names = "mu.x_10"))),
               "Locations of observations without measurement error information")
  fit <- fit_quiet(fd, spatialCovs = s$covs, par = p, prior = data.frame(mean = fd$x[11], sd = 1, row.names = "mu.x_11"))
  expect_s3_class(fit, "fitLangevin")
  # a base-name prior on 'mu' only applies to the estimated locations
  fit2 <- fit_quiet(fd, spatialCovs = s$covs, par = p, prior = data.frame(mean = 0, sd = 1e4, row.names = "mu"))
  expect_equal(fit2$estimates$random$mu$est$mu.x[10], fd$x[10])
})

test_that("a user map on mu is combined with the known locations", {
  s <- known_sim(n = 30, seed = 8)
  d <- s$dat; d$x.err <- d$y.err <- 2; d$x.err[10] <- d$y.err[10] <- NA
  fd <- formatData(d, time.unit = "hours")
  n <- nrow(fd); um <- seq_len(2 * n); um[c(29, 30)] <- NA                    # user fixes state 15 (at par$mu)
  pm <- cbind(fd$x, fd$y)
  fit <- fit_quiet(fd, spatialCovs = s$covs, par = list(beta = 0, sigma = 5, gamma = 0.5, mu = pm, vel = matrix(0, n, 2)),
                   map = list(mu = factor(um)))
  m <- matrix(as.character(fit$tmb_setup$map$mu), nrow = 2)
  expect_equal(which(is.na(m[1, ])), c(10, 15))
  expect_equal(nlevels(fit$tmb_setup$map$mu), 2 * n - 4)
  mu <- fit$estimates$random$mu$est
  expect_equal(mu$mu.x[c(10, 15)], fd$x[c(10, 15)])
})

test_that("predLangevin, residuals and cAIC work with mixed known locations", {
  s <- known_sim(n = 40, seed = 9)
  d <- s$dat; d$x.err <- d$y.err <- 2; k <- c(6, 25)
  d$x.err[k] <- d$y.err[k] <- NA
  fd <- formatData(d, time.unit = "hours")
  fit <- fit_quiet(fd, spatialCovs = s$covs, par = list(beta = 0, sigma = 5, gamma = 0.5), getJointPrecision = TRUE)
  res <- suppressMessages(suppressWarnings(residuals(fit, fd, s$covs)))
  expect_equal(nrow(res), nrow(fd))
  expect_true(all(is.na(res$residual.x[k])))
  expect_true(all(is.finite(res$residual.x[-c(1, k)])))
  ca <- suppressMessages(suppressWarnings(cAIC(fit, fd, s$covs)))
  expect_true(all(is.finite(unlist(ca[vapply(ca, is.numeric, TRUE)]))))
  # prediction: pad the data with prediction times (formatData(predTimes=)); the fit's known locations stay fixed
  fdp <- formatData(d, predTimes = data.frame(id = "A", date = fd$date[6] + c(-1800, 1800)), time.unit = "hours")
  pr <- suppressMessages(suppressWarnings(predLangevin(fit, fdp, s$covs)))
  pm <- pr$estimates$random$mu$est
  kp <- which(is.na(fdp$x.err) & !is.na(fdp$x))
  expect_equal(pm$mu.x[kp], fdp$x[kp])
  expect_true(all(is.finite(pm$mu.x[is.na(fdp$x)])))
})
