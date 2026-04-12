#' Constructor and validator for dataLangevin class
#'
#' @param x A data frame containing the Langevin data.
#' @param time.unit Optional character string specifying the time unit.
#' @return A validated \code{dataLangevin} object.
#' @noRd
class_dataLangevin <- function(x, time.unit = NULL) {
  if (!is.data.frame(x)) stop("Object must be a data frame to be a 'dataLangevin' object.")

  # enforce universal required columns
  req_cols <- c("id", "date", "dt", "x", "y", "smaj", "smin", "eor", "x.err", "y.err")
  missing_cols <- setdiff(req_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for 'dataLangevin': ", paste(missing_cols, collapse = ", "))
  }

  # validate 'date' flexibility (POSIXt/Date from formatData, or numeric from simLangevin)
  if (!(inherits(x$date, "POSIXt") || inherits(x$date, "Date") || is.numeric(x$date))) {
    stop("The 'date' column must be of class POSIXt, Date, or numeric.")
  }

  # apply time.unit attribute if provided
  if (!is.null(time.unit)) {
    attr(x, "time.unit") <- time.unit
  }

  class(x) <- unique(c("dataLangevin", class(x)))

  return(x)
}
