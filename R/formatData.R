#' Format data for \code{\link{fitLangevin}}
#'
#' \code{formatData} takes a data frame or \code{sf} object containing tracking data and formats it for use with the \code{\link{fitLangevin}} function. It is a modified version of the \code{format_data} function from the \href{https://ianjonsen.github.io/aniMotum/}{aniMotum} package.
#'
#' @param data A data frame or \code{\link[sf]{sf}} object containing the tracking data to be formatted. The data must contain columns for the animal ID (``id''), date/time (``date''), location quality class (``lc''; if applicable), coordinates (``coord''), and measurement error (``epar'' and/or ``sderr''; if applicable). The column names can be specified using the corresponding arguments of this function.
#' @param id Character string specifying the column name for the animal ID. Default: ``id''.
#' @param date Character string specifying the column name for the date/time, which must be of class \code{\link{DateTimeClasses}}. Default: ``date''.
#' @param coord Character vector of length 2 specifying the column names for the (projected) coordinates. Default: c("x", "y").
#' @param lc Character string specifying the column name for the location quality class. Default: ``lc''.
#' @param epar Character vector of length 3 specifying the column names for the error ellipse parameters (semi-major axis, semi-minor axis, and error ellipse orientation). See Details. Default: c("smaj", "smin", "eor").
#' @param sderr Character vector of length 2 specifying the column names for the standard deviations of the error for the x and y coordinates. Default: c("x.sd", "y.sd").
#' @param emf An optional data frame containing error multiplication factors (or standard deviations) for location quality classes. Must contain columns \code{lc}, \code{emf.x}, and \code{emf.y} (see \code{\link{get_emf}}). If provided, these values will be used to fill in \code{x.sd} and \code{y.sd} for observations where neither error ellipse (``epar'') nor standard deviations (``sderr'') information is provided in \code{data}. Default: \code{NULL}.
#' @param time.unit Character string specifying the time unit for the time steps. Default: ``hours''.
#' @param tz Character string specifying the time zone for the date/time column. Default: ``UTC''.
#' @return A data frame of class \code{dataLangevin} containing the formatted tracking data. The data frame contains the following columns:
#' \item{id}{Animal ID}
#' \item{date}{Date/time of observation}
#' \item{dt}{Numeric time step between observations (in \code{time.unit})}
#' \item{x}{x-coordinate of the location}
#' \item{y}{y-coordinate of the location}
#' \item{lc}{Location quality class (if applicable)}
#' \item{smaj}{Semi-major axis of the error ellipse}
#' \item{smin}{Semi-minor axis of the error ellipse}
#' \item{eor}{Error ellipse orientation (in radians)}
#' \item{x.sd}{Standard deviation of the x-coordinate error}
#' \item{y.sd}{Standard deviation of the y-coordinate error}
#'
#' @details
#' \strong{Coordinate Systems and Units:}
#' The coordinates \strong{must be projected to a metric coordinate system} (i.e., not longitude/latitude). For error ellipse data, the columns specified in \code{epar} should contain the semi-major axis (on the same scale as the coordinates), semi-minor axis (on the same scale as the coordinates), and error ellipse orientation (in \strong{degrees from north}). \code{formatData} automatically converts the orientation from degrees to radians.
#' For x- and y-axis error data, the columns specified in \code{sderr} should contain the standard deviations of the errors for the x and y coordinates (on the same scale as the coordinates).
#'
#' \strong{Location Classes and Measurement Errors:}
#' \code{formatData} is highly flexible in how it handles measurement error information. Observations can have error ellipses (``epar''), standard deviations (``sderr''), missing errors that are automatically filled based on location quality class (``lc'') using \code{emf} (see \code{\link{get_emf}}), or no measurement error. Location quality classes (``lc'') are used to determine which observations have unknown errors that need to be filled using \code{emf}. These classes can include Argos Least Squares (3, 2, 1, 0, A, B, Z), GPS or Argos Kalman Filter (G), and Generic Locations (GL). The handling of measurement error for each type of location quality class is as follows:
#' \itemize{
#'   \item \strong{Argos Least Squares (3, 2, 1, 0, A, B, Z):} Typically lack explicit error parameters. If missing, these can be filled using the \code{emf} argument (see \code{\link{get_emf}}). However, if users explicitly provide error ellipses or standard deviations for any given observation in \code{data}, the user-specified values are preserved and \code{emf} filling is bypassed for this particular observation.
#'   \item \strong{GPS or Argos Kalman Filter (G):} Can be provided with error ellipses, standard deviations, or no errors (which will be filled by the \code{emf} table if provided).
#'   \item \strong{Generic Locations (GL):} Locations where standard deviations (\code{x.sd}, \code{y.sd}) are explicitly provided by the user. If \code{emf} is specified but no \code{x.sd} and \code{y.sd} are provided, these rows are not filled unless \code{emf} includes an entry for \code{GL} (otherwise an error is returned).
#'   \item \strong{No Measurement Error:} If the \code{emf} argument is \code{NULL} (the default), observations with \code{NA} for all \code{epar} columns and \code{NA} for all \code{sderr} columns are assumed to have no measurement error. Their existing location class is preserved.
#' }
#'
#' @examples
#' # exampleDat included in package; see ?exampleDat for details
#' head(exampleDat)
#'
#' formatDat <- formatData(exampleDat, time.unit = "mins")
#'
#' \dontrun{
#' # exampleCovs included in package; see ?exampleCovs for details
#' fit <- fitLangevin(formatDat, spatialCovs = exampleCovs)
#' }
#' @importFrom dplyr rename arrange select all_of everything
#' @importFrom sf st_coordinates st_drop_geometry st_is_longlat
#' @export
formatData <- function(data, id = "id", date = "date", coord = c("x", "y"), lc = "lc", epar = c("smaj", "smin", "eor"), sderr = c("x.sd", "y.sd"), emf = NULL, time.unit = "hours", tz = "UTC"){

  if (inherits(data, "sf")) {
    if (sf::st_is_longlat(data)) {
      stop("The provided 'sf' object uses unprojected longitude/latitude coordinates. Please project your data to a metric coordinate system (e.g., using sf::st_transform) before formatting.")
    }

    coords_mat <- sf::st_coordinates(data)
    data <- sf::st_drop_geometry(data)
    data[[coord[1]]] <- coords_mat[, 1]
    data[[coord[2]]] <- coords_mat[, 2]

  } else {
    x_vals <- data[[coord[1]]]
    y_vals <- data[[coord[2]]]

    if (all(x_vals >= -180 & x_vals <= 360, na.rm = TRUE) &&
        all(y_vals >= -90 & y_vals <= 90, na.rm = TRUE)) {
      warning("Coordinates appear to be unprojected longitude/latitude. Langevin models require projected coordinates (e.g., UTM). Please ensure your data is projected.")
    }
  }

  if(!inherits(data[[date]],"POSIXt") && !is.character(data[[date]])) {
    stop("data$date must be of class 'POSIXt' or a character string parseable to a date.")
  }

  if (!(lc %in% names(data))) {
    data[[lc]] <- NA_character_
  }

  if (!all(is.na(data[[lc]]))) {
    valid_lcs <- c("3","2","1","0","A","B","Z","G","GL")
    provided_lcs <- data[[lc]][!is.na(data[[lc]])]
    if (!all(provided_lcs %in% valid_lcs)) {
      stop("Invalid location classes detected. Allowed values: 3, 2, 1, 0, A, B, Z, G, GL.")
    }
  }

  out <- format_data(x = data, id = id, date = date, coord = coord, lc = lc, epar = epar, sderr = sderr, tz = tz)

  out <- out %>%
    dplyr::rename(x = dplyr::all_of(coord[1]), y = dplyr::all_of(coord[2])) %>%
    dplyr::arrange(id, date)

  dt_list <- lapply(split(out$date, out$id), function(t) {
    c(0, as.numeric(difftime(t[-1], t[-length(t)], units = time.unit)))
  })
  out$dt <- unsplit(dt_list, out$id)

  out <- out %>% dplyr::select(id, date, dt, x, y, lc, smaj, smin, eor, x.sd, y.sd, dplyr::everything())

  # --- EMF INTEGRATION ---
  if (!is.null(emf)) {
    if (!is.data.frame(emf) || !all(c("lc", "emf.x", "emf.y") %in% names(emf))) {
      stop("'emf' must be a data frame containing 'lc', 'emf.x', and 'emf.y' columns.")
    }

    # Identify rows that have valid coordinates but lack ANY error information
    missing_err_idx <- which(!is.na(out$x) & !is.na(out$y) &
                               is.na(out$smaj) & is.na(out$smin) & is.na(out$eor) &
                               is.na(out$x.sd) & is.na(out$y.sd))

    if (length(missing_err_idx) > 0) {

      # Determine which LCs actually need filling
      lcs_to_fill <- unique(as.character(out$lc[missing_err_idx]))
      missing_lcs <- setdiff(lcs_to_fill, as.character(emf$lc))

      # Only stop if an LC that NEEDS filling is missing from the EMF table
      if (length(missing_lcs) > 0) {
        stop("The following location classes require error filling but are missing from the 'emf' table: ", paste(missing_lcs, collapse = ", "))
      }

      match_idx <- match(as.character(out$lc[missing_err_idx]), as.character(emf$lc))
      valid_matches <- !is.na(match_idx)
      update_idx <- missing_err_idx[valid_matches]

      out$x.sd[update_idx] <- emf$emf.x[match_idx[valid_matches]]
      out$y.sd[update_idx] <- emf$emf.y[match_idx[valid_matches]]
    }
  }

  if(any(!is.na(out$eor)) && isFALSE(any(out$eor > pi, na.rm = TRUE))) {
    warning(epar[3], " values were converted to radians, but they appear to have been provided in radians rather than degrees from north. Please ensure that the ", epar[3], " column was provided in degrees from north.")
  }

  out$eor <- out$eor * pi / 180 # convert from degrees to radians

  out$id <- as.factor(out$id)
  out$lc <- as.factor(out$lc)
  rownames(out) <- NULL

  checkErrorData(out, coord = c("x", "y"))

  attr(out, "time.unit") <- time.unit
  if(!inherits(out, "dataLangevin")) class(out) <- append(class(out), "dataLangevin")

  return(out)
}
