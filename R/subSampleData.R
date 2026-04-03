#' Subsample tracking data
#'
#' This function subsamples tracking data by a specified sampling rate and randomly introduces missing values in specified columns. It ensures that the first observation of each track is always included in the subsample.
#'
#' @param data A \code{dataLangevin} object (as returned by \code{\link{simLangevin}} and \code{\link{formatData}}) containing tracking data to subsample.
#' @param samplingRate A numeric value specifying the desired sampling rate. For example, a value of 2 will on average keep every second observation. Must be \code{>=1}. Default is 1 (no subsampling).
#' @param propMissing A numeric value between 0 and 1 specifying the proportion of observations to randomly set as missing (NA) in the subsample. Default is 0 (no missing values).
#' @param col_to_na A character vector of column names in the data frame for which to introduce missing values. Default is \code{c("x", "y", "smaj", "smin", "eor", "x.sd", "y.sd")}. Ignored unless \code{propMissing > 0}.
#' @return A subsampled \code{dataLangevin} object with the same structure as \code{data}, but with fewer rows (if \code{samplingRate>1}) and missing values (if \code{propMissing>0}) in the \code{col_to_na} columns.
#' @examples
#' # subsample with a sampling rate of 10 and 10% missing values
#' # exampleDat is an example dataLangevin object included in the package
#' exampleSub <- subSampleData(exampleDat, samplingRate =  10, propMissing = 0.1)
#'
#' @importFrom stats rbinom
#' @export
subSampleData <- function(data, samplingRate = 1, propMissing = 0, col_to_na = c("x", "y", "smaj", "smin", "eor", "x.sd", "y.sd")){

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  if(!is.numeric(samplingRate) || samplingRate < 1) stop("samplingRate should be an integer >= 1")
  if(!is.numeric(propMissing) || propMissing < 0 || propMissing >= 1) stop("propMissing must be a numeric value >=0 and <1")

  cols_to_na <- col_to_na[col_to_na %in% names(data)]
  if(length(cols_to_na) == 0 && propMissing > 0) warning("No valid columns specified in 'col_to_na' for introducing missing values. No missing values will be introduced.")

  n_rows <- nrow(data)
  current_ids <- data$id

  first_idx <- which(!duplicated(current_ids)) # Always include the first observation of each track

  n_sample <- ceiling(n_rows / max(samplingRate, 1))
  n_remaining <- n_sample - length(first_idx)

  # Uniformly sample the rest
  pool_idx <- which(duplicated(current_ids))
  rand_idx <- sample(pool_idx, n_remaining, replace = FALSE)

  # Sort indices to maintain timeline order
  sampled_idx <- sort(c(first_idx, rand_idx))

  sub_dat <- data[sampled_idx, ]
  row.names(sub_dat) <- NULL

  time.unit <- attr(data, "time.unit")

  dt_list <- lapply(split(sub_dat$date, sub_dat$id), function(t) {
    if (!is.null(time.unit) && (inherits(t, "POSIXt") || inherits(t, "Date"))) {
      c(0, as.numeric(difftime(t[-1], t[-length(t)], units = time.unit)))
    } else {
      c(0, as.numeric(diff(t)))
    }
  })

  sub_dat$dt <- unsplit(dt_list, sub_dat$id)
  # -------------------------------------------------------------------

  sub_ids <- sub_dat$id
  pool_sub <- which(duplicated(sub_ids))

  n_missing <- stats::rbinom(n = 1, size = length(pool_sub), prob = propMissing)

  if(n_missing > 0) {
    na_idx <- sample(pool_sub, n_missing, replace = FALSE)

    for(col in cols_to_na) {
      sub_dat[[col]][na_idx] <- NA
    }
  }

  attr(sub_dat, "time.unit") <- time.unit
  class(sub_dat) <- class(data)

  return(sub_dat)
}
