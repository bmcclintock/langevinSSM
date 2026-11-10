# tests/testthat/test_ncores.R
# Multi-threading (OpenMP) via 'ncores'.  At most 2 threads are used (CRAN policy).  Without OpenMP support,
# ncores = 2 falls back to one thread, and the comparisons hold trivially.

test_that("ncores is validated", {
  expect_error(langevinSSM:::.check_ncores(0), "positive integer")
  expect_error(langevinSSM:::.check_ncores(1.5), "positive integer")
  expect_error(langevinSSM:::.check_ncores(c(1, 2)), "positive integer")
  expect_identical(langevinSSM:::.check_ncores(1), 1L)
  expect_identical(suppressMessages(langevinSSM:::.check_ncores(2, has_openmp = FALSE)), 1L)
  expect_true(langevinSSM:::.check_ncores(2) %in% 1:2)
})

test_that("has_omp reports the OpenMP support of both compiled parts", {
  h <- has_omp()
  expect_type(h, "logical")
  expect_length(h, 1L)
  expect_identical(as.logical(h), attr(h, "likelihood") && attr(h, "simulator"))
  expect_identical(attr(h, "likelihood"), langevinSSM:::.tmb_has_openmp())
  expect_identical(attr(h, "simulator"), langevinSSM:::.rcpp_has_openmp())
  if (attr(h, "simulator")) expect_true(attr(h, "nprocs") >= 1L) else expect_true(is.na(attr(h, "nprocs")))
  expect_message(has_omp(verbose = TRUE), "OpenMP support")
})

test_that("simLangevin gives identical tracks for any ncores", {
  par <- list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5)
  for (m in c("underdamped", "overdamped")) {
    p <- if (m == "overdamped") par[c("beta", "sigma")] else par
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
    s1 <- suppressMessages(simLangevin(m, spatialCovs = exampleCovs, par = p, nbAnimals = 3, obsPerAnimal = 200, ncores = 1))
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
    s2 <- suppressMessages(simLangevin(m, spatialCovs = exampleCovs, par = p, nbAnimals = 3, obsPerAnimal = 200, ncores = 2))
    expect_identical(s1, s2)
  }
})

test_that("fitLangevin, predLangevin and simLangevin.fitLangevin agree for ncores = 1 and 2", {
  skip_on_cran()
  par <- list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  sim <- suppressMessages(simLangevin(spatialCovs = exampleCovs, par = par, nbAnimals = 2, obsPerAnimal = 150,
                                      timeStep = 0.05,
                                      measurementError = list(smaj.sd = 1.5, smin.sd = 0.75, eor.lim = c(0, 180))))
  f1 <- suppressMessages(fitLangevin(sim, spatialCovs = exampleCovs, ncores = 1))
  f2 <- suppressMessages(fitLangevin(sim, spatialCovs = exampleCovs, ncores = 2))
  expect_equal(f1$objective, f2$objective, tolerance = 1e-8)
  expect_equal(f1$par, f2$par, tolerance = 1e-5)
  expect_equal(f1$estimates$natural[, "Std. Error"], f2$estimates$natural[, "Std. Error"], tolerance = 1e-4)
  # the number of threads of the TMB library is restored on exit
  expect_equal(as.integer(TMB::openmp(DLL = "langevinSSM_TMBExports")), 1L)

  p1 <- suppressMessages(predLangevin(f1, data = sim, spatialCovs = exampleCovs, calcSE = FALSE, ncores = 1))
  p2 <- suppressMessages(predLangevin(f1, data = sim, spatialCovs = exampleCovs, calcSE = FALSE, ncores = 2))
  expect_equal(p1$objective, p2$objective, tolerance = 1e-8)

  set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
  s1 <- suppressMessages(simLangevin(f1, data = sim, spatialCovs = exampleCovs, ncores = 1))
  set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
  s2 <- suppressMessages(simLangevin(f1, data = sim, spatialCovs = exampleCovs, ncores = 2))
  expect_equal(s1, s2, tolerance = 1e-6)
  expect_equal(nlevels(s1$id), 2L)
})
