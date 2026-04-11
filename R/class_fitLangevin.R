#' Constructor and validator for fitLangevin class
#'
#' @param x A list containing the fitted model output.
#' @return A validated \code{fitLangevin} object.
#' @noRd
class_fitLangevin <- function(x) {
  if (!is.list(x)) stop("x must be a list to be a 'fitLangevin' object.")

  # enforce required TMB/nlminb output elements
  req_elements <- c("par", "objective", "convergence", "message",
                    "iterations", "evaluations", "elapsedTime",
                    "estimates", "covariance", "conditions")
  missing_elements <- setdiff(req_elements, names(x))
  if (length(missing_elements) > 0) {
    stop("Missing required list elements for 'fitLangevin': ", paste(missing_elements, collapse = ", "))
  }

  class(x) <- unique(c("fitLangevin", class(x)))

  return(x)
}
