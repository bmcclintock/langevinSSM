# --- Helper: Constructor/validator for hierLangevin class ---
class_hierLangevin <- function(x) {
  if (!is.list(x)) stop("x must be a list to be a 'hierLangevin' object.")
  req_elements <- c("par", "objective", "convergence", "estimates", "covariance", "conditions", "signatures")
  missing_elements <- setdiff(req_elements, names(x))
  if (length(missing_elements) > 0) {
    stop("Missing required list elements for 'hierLangevin': ", paste(missing_elements, collapse = ", "))
  }
  class(x) <- unique(c("hierLangevin", class(x)))
  return(x)
}

#' Stage II of the two-stage hierarchical fitting approach
#'
#' @description
#' Combines the per-individual point estimates and covariance matrices from a
#' list of separate \code{\link{fitLangevin}} models into a multivariate linear
#' mixed model, fit via TMB. This follows Stage II of the multistage
#' Conditionally Independent Hierarchical Model (CIHM) approach of Johnson,
#' Brost & Hooten (2022, \emph{JABES} 27:382-400). It yields population-level
#' estimates of the habitat-selection coefficients (\code{mu_beta}), movement
#' parameters (\code{mu_sigma}, and for the underdamped model, \code{mu_gamma}),
#' and their corresponding between-individual standard deviations (\code{sd_beta},
#' \code{sd_sigma}, \code{sd_gamma}).
#'
#' @details
#' Writing \eqn{\hat\theta_i} for individual \eqn{i}'s Stage I working-scale
#' point estimate and \eqn{\hat S_i} for its Hessian-based covariance matrix,
#' Stage II fits the model:
#' \deqn{\hat\theta_i \sim \mathrm{MVN}(\mu + u_i, \hat S_i), \quad u_i \sim \mathrm{N}(0, D)}
#' with \eqn{\hat S_i} treated as \strong{fixed/known} and \eqn{D} a diagonal
#' matrix of between-individual variances.
#'
#' \strong{Maximum Penalized Likelihood (MPL):} Optional penalty
#' (Chung, Rabe-Hesketh, Dorie, Gelman & Liu 2013, Psychometrika) can be
#' applied to the population-level log standard deviations via the \code{mpl} argument to
#' help avoid boundary (zero-variance) collapse when a covariate's
#' between-individual signal is weak.
#'
#' @param fit_list A list of \code{fitLangevin} objects (one per individual).
#' Elements that are not valid \code{fitLangevin} objects (e.g., failed fits
#' returning \code{NULL} or \code{try-error}) are automatically dropped with a warning.
#' If the list is named, the names are used as individual IDs; otherwise, generic
#' IDs ("ind_1", "ind_2", etc.) are generated.
#' @param mpl Logical or numeric. Controls an optional maximum penalized likelihood
#' (MPL) penalty on the population-level log standard deviations
#' (Chung et al. 2013). The term added to the log-likelihood is
#' \code{(alpha - 1) * log(sd_j)} for each parameter \code{j}.
#' \itemize{
#'   \item \code{FALSE} (default): no penalty. Unbiased MLE, but risks boundary collapse.
#'   \item \code{1}: numerically identical to \code{FALSE} (zero penalty).
#'   \item \code{TRUE}: uses Chung et al.'s recommended default \code{alpha = 2}.
#'   \item Numeric scalar or vector of length \code{p}: sets \code{alpha} directly. Values \code{> 1} counteract boundary collapse.
#' }
#' @param id_col_name Character string. The column name to use for the individual IDs in the
#' returned \code{estimates$random} data frames. Default: \code{"id"}. Note this controls the
#' \emph{name} of the column. The \emph{values} within this column are populated by the names of
#' \code{fit_list} (if provided), or generated automatically (e.g., "ind_1", "ind_2").
#' @param silent Logical indicating whether or not to disable TMB tracing
#' information. Default: \code{FALSE}.
#' @param control List controlling the outer \code{\link[stats]{nlminb}}
#' optimization. Default: \code{list(trace = 0, iter.max = 1000, eval.max = 1000)}.
#'
#' @return An object of class \code{hierLangevin}, a list containing:
#' \item{par}{The optimized outer (fixed-effect) parameter vector.}
#' \item{objective}{The negative log-likelihood at the optimum.}
#' \item{convergence}{Convergence code from \code{\link[stats]{nlminb}} (0 = converged).}
#' \item{estimates}{A list with elements \code{working} (population-level
#' means and SDs on the working/log scale), \code{natural} (population-level
#' means on the natural scale, i.e. \code{mu_sigma}, \code{mu_gamma}
#' exponentiated, with delta-method SEs), and \code{random} (per-individual
#' natural-scale predictions/BLUPs and their SEs).}
#' \item{covariance}{A list with elements \code{working} and \code{natural} containing the Hessian-based covariance matrices for the population-level parameters. It also includes the sparse \code{jointPrecision} matrix of the random effects.}
#' \item{conditions}{A list recording the call's arguments for reproducibility.}
#' \item{signatures}{A list containing the spatial covariate signature shared by all individuals in the model.}
#'
#' @references
#' Johnson DS, Brost BM, Hooten MB (2022). Greater Than the Sum of its Parts:
#' Computationally Flexible Bayesian Hierarchical Modeling. \emph{Journal of
#' Agricultural, Biological and Environmental Statistics}, 27, 382-400.
#'
#' Lele SR, Glen K, Ponciano JM (2025). Practical Consequences of the Bias in
#' the Laplace Approximation to Marginal Likelihood for Hierarchical Models.
#' \emph{Entropy}, 27(3), 289.
#'
#' Chung Y, Rabe-Hesketh S, Dorie V, Gelman A, Liu J (2013). A nondegenerate
#' penalized likelihood estimator for variance parameters in multilevel
#' models. \emph{Psychometrika}, 78(4), 685-709.
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb sd
#' @export
hierLangevin <- function(fit_list, mpl = FALSE, id_col_name = "id", silent = FALSE,
                         control = list(trace = 0, iter.max = 1000, eval.max = 1000)) {

  if (!is.list(fit_list) || length(fit_list) == 0 || inherits(fit_list, "fitLangevin")) {
    stop("'fit_list' must be a list containing fitted 'fitLangevin' objects.")
  }

  # Filter out failed fits (e.g., try-errors or NULLs from lapply loops)
  valid_idx <- vapply(fit_list, inherits, logical(1), "fitLangevin")
  if (!all(valid_idx)) {
    warning(sprintf("Dropped %d element(s) from 'fit_list' that were not valid 'fitLangevin' objects.",
                    sum(!valid_idx)))
    fit_list <- fit_list[valid_idx]
  }

  n <- length(fit_list)
  if (n < 2) stop("hierLangevin requires at least 2 successful individual fits.")

  # Assign IDs
  ids <- names(fit_list)
  if (is.null(ids)) {
    ids <- paste0("ind_", seq_len(n))
  }

  # Extract parameter names and configuration from the first fit
  param_names <- rownames(fit_list[[1]]$estimates$working)
  p <- length(param_names)
  model_type <- fit_list[[1]]$conditions$model # Assuming fitLangevin stores the model type

  # Extract barrier and scaling metadata to pass downstream
  base_barrier <- fit_list[[1]]$conditions$barrier
  base_lambda <- fit_list[[1]]$conditions$lambda
  base_sf <- if (!is.null(fit_list[[1]]$conditions$scaleFactor)) fit_list[[1]]$conditions$scaleFactor else 1
  base_covs_sig <- fit_list[[1]]$signatures$covs

  # Initialize structures
  theta_hat <- matrix(NA_real_, n, p, dimnames = list(ids, param_names))
  Shat <- array(NA_real_, dim = c(p, p, n), dimnames = list(param_names, param_names, ids))

  # Populate Stage I estimates and covariance matrices
  for (i in seq_len(n)) {
    est_i <- fit_list[[i]]$estimates$working
    if (!identical(rownames(est_i), param_names)) {
      stop(sprintf("Working-scale parameter names for individual '%s' do not match individual '%s'. ",
                   ids[i], ids[1]),
           "All individuals must be fitted with the exact same model and spatial covariates.")
    }

    # Check for identical structural arguments
    if (!identical(fit_list[[i]]$conditions$barrier, base_barrier)) {
      stop(sprintf("Barrier argument for individual '%s' does not match individual '%s'. All individuals must be fitted with the same barrier constraint.", ids[i], ids[1]))
    }
    if (!identical(fit_list[[i]]$conditions$lambda, base_lambda)) {
      stop(sprintf("Lambda penalty for individual '%s' does not match individual '%s'. All individuals must be fitted with the same lambda.", ids[i], ids[1]))
    }
    if (!identical(fit_list[[i]]$signatures$covs, base_covs_sig)) {
      if (!isTRUE(all.equal(fit_list[[i]]$signatures$covs, base_covs_sig, tolerance = 1e-5))) {
        stop(sprintf("spatialCovs signature for individual '%s' does not match individual '%s'. All individuals must be fitted with the exact same spatial covariates.", ids[i], ids[1]))
      }
    }
    sf_i <- if (!is.null(fit_list[[i]]$conditions$scaleFactor)) fit_list[[i]]$conditions$scaleFactor else 1
    if (!identical(sf_i, base_sf)) {
      stop(sprintf("scaleFactor argument for individual '%s' does not match individual '%s'. All individuals must be fitted with the same scaleFactor.", ids[i], ids[1]))
    }

    theta_hat[i, ] <- est_i$Estimate

    cov_i <- fit_list[[i]]$covariance$working
    if (is.null(cov_i) || any(dim(cov_i) != c(p, p))) {
      stop(sprintf("Working-scale covariance matrix unavailable for individual '%s'. ", ids[i]),
           "Ensure calcSE = TRUE was used for every individual fit in Stage I.")
    }
    Shat[, , i] <- as.matrix(cov_i)
  }

  # MPL alpha handling
  if (is.numeric(mpl) && length(mpl) >= 1 && any(mpl < 1, na.rm = TRUE)) {
    warning("'mpl' alpha < 1 applies a penalty of (alpha - 1) * log(sd), which pushes the ",
            "between-individual standard deviations TOWARDS zero and makes boundary collapse ",
            "more likely. Use mpl = 1 (or FALSE) for no penalty, and mpl > 1 to counteract ",
            "boundary collapse.", call. = FALSE)
  }
  if (isFALSE(mpl)) {
    use_mpl <- 0L
    mpl_alpha <- rep(2, p)
  } else if (isTRUE(mpl)) {
    use_mpl <- 1L
    mpl_alpha <- rep(2, p)
  } else if (is.numeric(mpl)) {
    use_mpl <- 1L
    if (length(mpl) == 1) {
      mpl_alpha <- rep(mpl[1], p)
    } else if (length(mpl) == p) {
      mpl_alpha <- as.numeric(mpl)
    } else {
      stop(sprintf("'mpl' must be a scalar or a vector of length %d (number of working-scale parameters), got length %d", p, length(mpl)))
    }
  } else {
    stop("'mpl' must be FALSE, TRUE, or a numeric scalar/vector of alpha values.")
  }

  # Determine which columns are log-scale movement parameters (sigma/gamma)
  col_type <- as.integer(param_names %in% c("log_sigma", "log_gamma"))
  col_shift <- rep(0, p)

  # data_list matches the TMB template exactly
  data_list <- list(theta_hat = theta_hat, Shat = Shat,
                    use_mpl = use_mpl, mpl_alpha = mpl_alpha,
                    col_type = col_type, col_shift = col_shift)

  par_list <- list(
    mu = as.numeric(colMeans(theta_hat)),
    log_sd = log(pmax(apply(theta_hat, 2, sd), 1e-3)),
    u = matrix(0, n, p)
  )

  obj <- TMB::MakeADFun(c(model = "hierLangevin", data_list), par_list,
                        random = "u", DLL = "langevinSSM_TMBExports", silent = silent)

  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = control)

  if (opt$convergence != 0) {
    warning(sprintf("hierLangevin outer optimization failed to converge (nlminb code: %s, message: %s)",
                    opt$convergence, opt$message))
  }

  # Ensure getJointPrecision = TRUE is called here so we don't need to rebuild tapes later
  sdr <- try(TMB::sdreport(obj, getJointPrecision = TRUE), silent = TRUE)

  estimates <- list()
  covariance <- list()

  if (inherits(sdr, "try-error")) {
    warning("TMB::sdreport failed to calculate standard errors for hierLangevin. Point estimates are available, but SEs are unavailable.")
    obj$fn(opt$par)
    rep_vals <- obj$report()
    working_names <- c(paste0("mu_", param_names), paste0("sd_", param_names))
    estimates$working <- data.frame(Estimate = c(rep_vals$mu, rep_vals$sd),
                                    "Std. Error" = NA_real_, check.names = FALSE,
                                    row.names = working_names)
    nat_names <- paste0("mu_", ifelse(col_type == 1, sub("^log_", "", param_names), param_names))
    estimates$natural <- data.frame(Estimate = rep_vals$mu_nat, "Std. Error" = NA_real_,
                                    check.names = FALSE, row.names = nat_names)
    ind_names <- ifelse(col_type == 1, sub("^log_", "", param_names), param_names)
    estimates$random <- list()
    for (j in seq_len(p)) {
      estimates$random[[ind_names[j]]] <- list(
        est = data.frame(id = ids, est = rep_vals$ind_nat[, j]),
        se = NULL
      )
      names(estimates$random[[ind_names[j]]]$est)[1] <- id_col_name
    }
  } else {
    rep_summary <- summary(sdr, "report")

    mu_rows <- rep_summary[rownames(rep_summary) == "mu", , drop = FALSE]
    sd_rows <- rep_summary[rownames(rep_summary) == "sd", , drop = FALSE]
    rownames(mu_rows) <- paste0("mu_", param_names)
    rownames(sd_rows) <- paste0("sd_", param_names)
    estimates$working <- as.data.frame(rbind(mu_rows, sd_rows))

    mu_nat_rows <- rep_summary[rownames(rep_summary) == "mu_nat", , drop = FALSE]
    nat_names <- paste0("mu_", ifelse(col_type == 1, sub("^log_", "", param_names), param_names))
    rownames(mu_nat_rows) <- nat_names
    estimates$natural <- as.data.frame(mu_nat_rows)

    ind_nat_rows <- rep_summary[rownames(rep_summary) == "ind_nat", , drop = FALSE]
    ind_names <- ifelse(col_type == 1, sub("^log_", "", param_names), param_names)
    estimates$random <- list()

    # NOTE: TMB's ADREPORT of a matrix flattens in column-major order
    for (j in seq_len(p)) {
      idx <- ((j - 1) * n + 1):(j * n)
      estimates$random[[ind_names[j]]] <- list(
        est = data.frame(id = ids, est = ind_nat_rows[idx, "Estimate"]),
        se  = data.frame(id = ids, se  = ind_nat_rows[idx, "Std. Error"])
      )
      names(estimates$random[[ind_names[j]]]$est)[1] <- id_col_name
      names(estimates$random[[ind_names[j]]]$se)[1]  <- id_col_name
    }

    # Store covariance matrices
    covariance$natural <- sdr$cov
    covariance$working <- sdr$cov.fixed
    if (!is.null(sdr$jointPrecision)) {
      covariance$random <- list(jointPrecision = sdr$jointPrecision)
    }

    # Ensure working covariance matrix carries exact dimension names matching estimates
    working_names <- c(paste0("mu_", param_names), paste0("sd_", param_names))
    rownames(covariance$working) <- working_names
    colnames(covariance$working) <- working_names
  }

  tmb_setup <- list(
    parList = obj$env$parList(opt$par)
  )

  out <- list(
    par = opt$par,
    objective = opt$objective,
    convergence = opt$convergence,
    message = opt$message,
    iterations = opt$iterations,
    evaluations = opt$evaluations,
    estimates = estimates,
    covariance = covariance,
    conditions = list(mpl = mpl, silent = silent, control = control,
                      model = model_type, param_names = param_names, ids = ids,
                      barrier = base_barrier, lambda = base_lambda, scaleFactor = base_sf),
    signatures = list(covs = base_covs_sig),
    tmb_setup = tmb_setup,
    sdreport = if (!inherits(sdr, "try-error")) sdr else NULL
  )

  class_hierLangevin(out)
}
