# tests/testthat/test_gof_tests.R
# Unit tests of the goodness-of-fit tests attached by residuals.fitLangevin() (no model fit needed).

mk_res <- function(rx, ry, id) {
  structure(data.frame(id = id, date = seq_along(id), residual.x = rx, residual.y = ry),
            class = c("resLangevin", "data.frame"))
}
gof <- function(...) suppressMessages(langevinSSM:::gof_tests(...))

test_that("gof_tests returns the reduced set of tests with effect sizes", {
  set.seed(1)
  n <- 600
  tt <- gof(mk_res(rnorm(n), rnorm(n), rep(1:3, each = n / 3)), nsim = 200)
  expect_equal(tt$metric, c("KS_x", "KS_y", "KS_mah_shape", "Disp_mah", "LB_x", "LB_y", "LB_mah"))
  expect_named(tt, c("metric", "statistic", "p.value", "effect"))
  expect_true(all(tt$p.value >= 0 & tt$p.value <= 1))
  expect_false("KS_mah" %in% tt$metric)
  # effect = KS distance for the KS tests and variance error for Disp_mah
  expect_equal(tt$effect[1:3], tt$statistic[1:3])
  expect_equal(tt$effect[4], tt$statistic[4] - 1)
})

test_that("Ljung-Box tests are computed within tracks and pooled exactly", {
  set.seed(2)
  id <- rep(c("a", "b", "c"), times = c(200, 350, 90))
  rx <- rnorm(length(id)); ry <- rnorm(length(id))
  tt <- gof(mk_res(rx, ry, id), nsim = 50)
  Q <- 0; df <- 0
  for (g in unique(id)) {
    z <- rx[id == g]; L <- floor(log(length(z)))
    b <- stats::Box.test(z, lag = L, type = "Ljung-Box")
    Q <- Q + b$statistic; df <- df + L
  }
  expect_equal(tt$statistic[tt$metric == "LB_x"], unname(Q), tolerance = 1e-10)
  expect_equal(tt$p.value[tt$metric == "LB_x"], unname(stats::pchisq(Q, df, lower.tail = FALSE)), tolerance = 1e-10)
})

test_that("Ljung-Box tests ignore the jump between consecutive tracks", {
  # each track is white noise around its own (very different) mean: concatenating the tracks
  # would create strong spurious autocorrelation, the per-track test must not see it
  set.seed(3)
  id <- rep(1:6, each = 150)
  rx <- rnorm(length(id)) + rep(c(-3, 3, -3, 3, -3, 3), each = 150)
  ry <- rnorm(length(id))
  tt <- gof(mk_res(rx, ry, id), nsim = 50)
  expect_gt(stats::Box.test(rx, lag = floor(log(length(rx))), type = "Ljung-Box")$statistic, 1000)
  expect_gt(tt$p.value[tt$metric == "LB_x"], 0.001)
  expect_lt(tt$effect[tt$metric == "LB_x"], 0.1)
})

test_that("Ljung-Box effect size reflects the autocorrelation", {
  set.seed(4)
  id <- rep(1:2, each = 3000)
  rx <- c(stats::arima.sim(list(ar = 0.3), 3000), stats::arima.sim(list(ar = 0.3), 3000))
  rx <- rx / stats::sd(rx)
  tt <- gof(mk_res(rx, rnorm(length(id)), id), nsim = 50)
  lbx <- tt[tt$metric == "LB_x", ]
  expect_lt(lbx$p.value, 1e-10)
  expect_equal(lbx$effect, 0.3, tolerance = 0.15)
  expect_lt(tt$effect[tt$metric == "LB_y"], 0.06)
})

test_that("gof_tests handles very short tracks and NA residuals", {
  set.seed(5)
  id <- c(rep(1, 2), rep(2, 300))
  rx <- rnorm(302); ry <- rnorm(302); rx[c(10, 20)] <- NA
  tt <- gof(mk_res(rx, ry, id), nsim = 50)
  expect_true(all(is.finite(tt$p.value)))
})
