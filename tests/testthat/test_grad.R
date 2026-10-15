test_that("grad option is validated, stored and passed to the template", {
  expect_identical(langevinSSM:::.grad_degree("bilinear"), 1L)
  expect_identical(langevinSSM:::.grad_degree("cubic"), 3L)
  expect_identical(langevinSSM:::.grad_degree("quintic"), 5L)
  expect_error(langevinSSM:::.grad_degree("spline"))
  expect_identical(langevinSSM:::.cond_grad(list()), "quintic")
  expect_identical(langevinSSM:::.cond_grad(list(grad = "cubic")), "cubic")
})

test_that("simLangevin and fitLangevin run with every grad option", {
  skip_on_cran()
  set.seed(1)
  ec <- langevinSSM::exampleCovs
  par <- list(beta = c(-4, 6, 5, -0.1), sigma = 0.2, gamma = 0.02)
  ip <- c(mean(terra::ext(ec[[1]])[1:2]), mean(terra::ext(ec[[1]])[3:4]))
  for (g in c("bilinear", "cubic", "quintic")) {
    set.seed(2)
    s <- suppressMessages(simLangevin(model = "underdamped", spatialCovs = ec, par = par, grad = g,
                                      timeStep = 0.5, obsPerAnimal = 60, initialPosition = ip,
                                      measurementError = list(smaj.sd = 0.5, smin.sd = 0.25)))
    expect_s3_class(s, "dataLangevin")
    f <- suppressWarnings(suppressMessages(fitLangevin(s, spatialCovs = ec, model = "underdamped", grad = g,
                                                       calcSE = FALSE, control = list(iter.max = 5, eval.max = 10))))
    expect_identical(f$conditions$grad, g)
    expect_true(is.finite(f$objective))
  }
  expect_error(suppressMessages(simLangevin(model = "underdamped", spatialCovs = ec, par = par, grad = "linear",
                                            obsPerAnimal = 10, initialPosition = ip)))
})
