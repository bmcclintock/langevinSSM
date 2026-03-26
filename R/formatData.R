#' Format data for \code{\link{fitLangevin}}
#'
#' \code{formatData} takes a data frame containing tracking data and formats it for use with the \code{\link{fitLangevin}} function. It is a modified version of the \code{format_data} function from the \href{https://ianjonsen.github.io/aniMotum/}{aniMotum} package.
#'
#' @param data A data frame containing the tracking data to be formatted. The data frame must contain columns for the animal ID (``id''), date/time (``date''), location quality class (``lc'', which can be 3, 2, 1, 0, A, B, or Z for Argos least squares and "G" for Argos error ellipse or GPS), coordinates, and error parameters (if applicable). The column names can be specified using the corresponding arguments of this function.
#' @param id Character string specifying the column name for the animal ID. Default: ``id''.
#' @param date Character string specifying the column name for the date/time, which must be of class \code{\link{DateTimeClasses}}. Default: ``date''.
#' @param coord Character vector of length 2 specifying the column names for the coordinates. Default: c("x", "y").
#' @param lc Character string specifying the column name for the location quality class. Default: ``lc''.
#' @param epar Character vector of length 3 specifying the column names for the error ellipse parameters (semi-major axis, semi-minor axis, and error ellipse orientation). Default: c("smaj", "smin", "eor").
#' @param sderr Character vector of length 2 specifying the column names for the standard deviations of the error for the x and y coordinates. Default: c("x.sd", "sd.y").
#' @param time.unit Character string specifying the time unit for the time steps. Default: ``hours''.
#' @param tz Character string specifying the time zone for the date/time column. Default: ``UTC''.
#' @return A data frame of class \code{dataLangevin} containing the formatted tracking data. The data frame contains the following columns:
#' \item{id}{Animal ID}
#' \item{date}{Date/time of observation}
#' \item{dt}{Time step between observations (in \code{time.unit})}
#' \item{x}{x-coordinate of the location}
#' \item{y}{y-coordinate of the location}
#' \item{lc}{Location quality class}
#' \item{smaj}{Semi-major axis of the error ellipse}
#' \item{smin}{Semi-minor axis of the error ellipse}
#' \item{eor}{Error ellipse orientation}
#' \item{x.sd}{Standard deviation of the x-coordinate error}
#' \item{y.sd}{Standard deviation of the y-coordinate error}
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
#' @export
formatData <- function(data, id = "id", date = "date", coord = c("x", "y"), lc = "lc", epar = c("smaj", "smin", "eor"), sderr = c("x.sd", "y.sd"), time.unit = "hours", tz = "UTC"){

  if(is.null(data[[lc]])) stop("data$lc is missing")
  if(!all(data[[lc]] %in% c(3,2,1,0,"A","B","Z","G"))) stop("data$lc can only be '3', '2', '1', '0', 'A', 'B', and 'Z' for location quality classes or 'G' for GPS data")
  if(!inherits(data[[date]],"POSIXt")) stop("data$date must be of class 'POSIXt'")
  out <- format_data(x = data, id = id, date = date, coord = coord, lc = lc, epar = epar, sderr = sderr, tz = tz)
  out <- out %>% dplyr::rename(x=coord[1],y=coord[2]) %>%
    dplyr::arrange(id,date) %>%
    dplyr::mutate(dt=do.call(c,mapply(function(x) c(0,diff(out$date[which(out$id==x)],units=time.unit)), unique(out$id), SIMPLIFY = FALSE))) %>%
    dplyr::select(id,date,dt,x,y,lc,smaj,smin,eor,x.sd,y.sd)
  out$id <- as.factor(out$id)
  out$lc <- as.factor(out$lc)
  rownames(out) <- NULL
  attr(out,"time.unit") <- time.unit
  if(!inherits(out,"dataLangevin")) class(out) <- append(class(out),"dataLangevin")
  return(out)
}
