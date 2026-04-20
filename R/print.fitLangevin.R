#' Print a fitLangevin object
#' @method print fitLangevin
#'
#' @param x A \code{fitLangevin} object returned by \code{\link{fitLangevin}}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @importFrom stats printCoefmat
#' @export
print.fitLangevin <- function(x, ...) {

  cat("\nHabitat-Driven Langevin Diffusion Model\n")
  cat("=======================================\n")

  # Determine model type
  model_type <- ifelse("vel" %in% names(x$estimates$random), "Underdamped", "Overdamped")
  cat("Model type:       ", model_type, "\n")

  # Convergence status
  conv_text <- ifelse(x$convergence == 0, "Successful", paste("Failed (Code", x$convergence, ")"))
  cat("Convergence:      ", conv_text, "\n")
  if (x$convergence != 0 && !is.null(x$message)) {
    cat("Message:          ", x$message, "\n")
  }

  cat("Max Log-Likelihood:", -x$objective, "\n")
  cat("Optimization time: ", round(x$elapsedTime[3], 2), "seconds\n\n")

  cat("Parameter Estimates (Natural Scale):\n")
  cat("---------------------------------------\n")

  # Clean up the rownames for the natural estimates matrix
  nat_est <- x$estimates$natural

  stats::printCoefmat(nat_est, digits = 4, signif.stars = FALSE, na.print = "NA", ...)

  boundsWarning(x)

  invisible(x)
}
