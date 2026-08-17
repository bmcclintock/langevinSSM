#' Calculate Conditional AIC for a fitted Langevin model
#'
#' Computes the Conditional Akaike Information Criterion (cAIC) using Method 2 from
#' Zheng, Cadigan, and Thorson (2024).
#'
#' @param fit A \code{fitLangevin} object.
#' @param data A \code{dataLangevin} object used to fit the model.
#' @param spatialCovs A list of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param nSims Integer. The number of Rademacher vectors used for Hutchinson's trace estimator (if the exact sparse inverse fails). Default: 200.
#'
#' @return A list containing the calculated \code{cAIC}, effective degrees of freedom (\code{EDF}),
#' the conditional negative log-likelihood (\code{NLL_cond}), the random effects trace penalty (\code{trace_penalty}),
#' the number of actively estimated fixed effects (\code{p}), and the total number of random effects (\code{q}).
#'
#' @references
#' Zheng, Y., Cadigan, N. G., & Thorson, J. T. (2024). A distributionally robust conditional
#' Akaike information criterion for random-effect models. arXiv preprint arXiv:2411.14185.
#'
#' @export
cAIC <- function(fit, data, spatialCovs, nSims = 200) {

  if (!requireNamespace("sparseinv", quietly = TRUE)) {
    stop("Package 'sparseinv' is required to calculate cAIC. Please install it.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required to calculate cAIC. Please install it.")
  }

  if (!inherits(fit, "fitLangevin")) {
    stop("'fit' must be a 'fitLangevin' object.")
  }

  if (is.null(fit$covariance$random$jointPrecision)) {
    stop("The model must be fitted with 'getJointPrecision = TRUE' to calculate cAIC. Please refit the model.")
  }

  verify_signatures(fit, data = data, spatialCovs = spatialCovs)

  scaleFactor <- fit$conditions$scaleFactor
  model <- fit$conditions$model
  lambda <- fit$conditions$lambda
  barrier <- fit$conditions$barrier
  coord <- fit$conditions$coord
  smoothGradient <- fit$conditions$smoothGradient
  npoints <- fit$conditions$npoints
  curweight <- fit$conditions$curweight
  zetaScale <- fit$conditions$zetaScale

  spatialCovs_copy <- spatialCovs
  barrier_sdf <- NULL
  if (!is.null(barrier)) {
    barrier_sdf <- spatialCovs_copy[[barrier]]
    spatialCovs_copy[[barrier]] <- NULL
  }

  dat_joint <- build_tmb_data(
    data = data,
    spatialCovs = spatialCovs_copy,
    model = model,
    coord = coord,
    scaleFactor = scaleFactor,
    smoothGradient = smoothGradient,
    npoints = npoints,
    curweight = curweight,
    zetaScale = zetaScale,
    barrier_sdf = barrier_sdf,
    lambda = lambda
  )

  dat_joint <- c(dat_joint, fit$tmb_setup$priors)

  parList_mle <- fit$tmb_setup$parList
  tmbmap <- fit$tmb_setup$map
  re <- fit$tmb_setup$random

  if (is.null(tmbmap)) tmbmap <- list()

  message("   Evaluating conditional log-likelihood...")

  obj_joint <- suppressWarnings(try({
    TMB::MakeADFun(
      data = c(model = "langevinSSM", dat_joint),
      par = parList_mle,
      map = tmbmap,
      random = re,
      DLL = "langevinSSM_TMBExports",
      silent = TRUE
    )
  }, silent = TRUE))

  if (inherits(obj_joint, "try-error")) {
    stop("Failed to construct the joint objective function: ", attr(obj_joint, "condition")$message)
  }

  obj_joint$fn(obj_joint$par)
  NLL_joint_states <- as.numeric(obj_joint$env$f(obj_joint$env$last.par))

  message("   Evaluating prior precision matrix...")

  dat_prior <- dat_joint
  dat_prior$isd[] <- 0

  obj_prior <- suppressWarnings(try({
    TMB::MakeADFun(
      data = c(model = "langevinSSM", dat_prior),
      par = parList_mle,
      map = tmbmap,
      random = re,
      DLL = "langevinSSM_TMBExports",
      silent = TRUE
    )
  }, silent = TRUE))

  if (inherits(obj_prior, "try-error")) {
    stop("Failed to construct the prior objective function: ", attr(obj_prior, "condition")$message)
  }

  obj_prior$fn(obj_prior$par)
  NLL_prior_states <- as.numeric(obj_prior$env$f(obj_prior$env$last.par))

  NLL_cond <- NLL_joint_states - NLL_prior_states

  sd_prior <- suppressWarnings(try({
    TMB::sdreport(obj_prior, getJointPrecision = TRUE)
  }, silent = TRUE))

  if (inherits(sd_prior, "try-error") || is.null(sd_prior$jointPrecision)) {
    stop("Failed to extract the prior joint precision matrix.")
  }

  Q_joint_full <- fit$covariance$random$jointPrecision
  Q_prior_full <- sd_prior$jointPrecision

  re_names <- c("mu", "vel")
  idx_joint <- which(rownames(Q_joint_full) %in% re_names)
  idx_prior <- which(rownames(Q_prior_full) %in% re_names)

  Q_joint_re <- Q_joint_full[idx_joint, idx_joint, drop = FALSE]
  Q_prior    <- Q_prior_full[idx_prior, idx_prior, drop = FALSE]

  q_re <- ncol(Q_prior)
  p_fe <- length(fit$par)

  if (q_re == 0) stop("The random effects precision matrix is empty (0x0).")
  if (ncol(Q_joint_re) != q_re) stop("Dimension mismatch between Q_prior and Q_joint_re.")

  message("   Calculating sparse inverse and trace penalty...")

  Q_joint_sym <- Matrix::forceSymmetric(Q_joint_re)

  trace_term <- tryCatch({
    invQ_sparse <- suppressWarnings(sparseinv::Takahashi_Davis(Q_joint_sym))
    sum(Q_prior * invQ_sparse)

  }, error = function(e) {
    message("      Exact sparse inverse subset failed. Executing Hutchinson's Trace Estimator...")
    tr_est <- 0
    chol_Q_joint <- Matrix::Cholesky(Q_joint_sym, super = TRUE)

    for (i in 1:nSims) {
      z <- matrix(sample(c(-1, 1), q_re, replace = TRUE), ncol = 1)
      y <- as.numeric(as.vector(Matrix::solve(chol_Q_joint, z)))
      Q_y <- as.numeric(as.vector(Q_prior %*% y))
      tr_est <- tr_est + sum(as.numeric(z) * Q_y)
    }
    return(tr_est / nSims)
  })

  EDF <- q_re - trace_term
  caic_val <- 2 * NLL_cond + 2 * (p_fe + EDF)

  res <- list(
    cAIC = caic_val,
    EDF = EDF,
    NLL_cond = NLL_cond,
    trace_penalty = trace_term,
    p = p_fe,
    q = q_re
  )

  class(res) <- "caicLangevin"
  return(res)
}

#' @export
print.caicLangevin <- function(x, ...) {
  cat("\nConditional Akaike Information Criterion (cAIC)\n")
  cat("===============================================\n")
  cat(sprintf("cAIC:                  %.2f\n", x$cAIC))
  cat(sprintf("Conditional NLL:       %.2f\n", x$NLL_cond))
  cat("-----------------------------------------------\n")
  cat(sprintf("Effective DF (EDF):    %.2f\n", x$EDF))
  cat(sprintf("Trace Penalty:         %.2f\n", x$trace_penalty))
  cat(sprintf("Fixed Effects (p):     %d\n", x$p))
  cat(sprintf("Random Effects (q):    %d\n\n", x$q))

  invisible(x)
}
