#' Plot One-Step-Ahead (OSA) Residuals for \code{fitLangevin} model
#'
#' Generates ggplot2 diagnostic plots (Q-Q plots and ACF plots) for the OSA residuals of a fitted Langevin diffusion model. Plots are made separately for the x and y observed location residuals (which should follow a normal distribution if the model fits well), as well as for the squared Mahalanobis distance of the observed location residuals (which should follow a Chi-Square distribution with 2 degrees of freedom if the model fits well).
#'
#' @param fit A \code{fitLangevin} object.
#' @param tracks Optional. Vector of track IDs to plot separately, or \code{"all"} to plot each track in the dataset individually. If \code{NULL} (default), residuals for all tracks are aggregated into a single set of plots.
#' @return List of \code{ggplot} objects (or a nested list of \code{ggplot} objects if plotting by track).
#' @examples
#' par <- list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5)
#' measurementError <- list(smaj.sd = 1.5, smin.sd = 0.75, eor = c(0,180))
#'
#' set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
#' # ridiculously small dataset for example purposes
#' smallDat <- simLangevin(par = par,
#'                           spatialCovs = exampleCovs,
#'                           nbAnimals = 3,
#'                           obsPerAnimal = 50,
#'                           measurementError = measurementError)
#'
#' fit <- fitLangevin(model = "underdamped",
#'                    data = smallDat,
#'                    spatialCovs = exampleCovs,
#'                    silent = TRUE,
#'                    control = list(trace = 1),
#'                    calcOSA = TRUE)
#'
#' plotResiduals(fit)
#'
#' plotResiduals(fit, tracks = c("1", "3"))
# #' @importFrom ggplot2 ggplot aes geom_hline geom_segment geom_abline labs theme_minimal geom_ribbon geom_point scale_color_manual theme
#' @importFrom stats acf qnorm qchisq ppoints quantile dnorm dchisq na.omit
#' @export
plotResiduals <- function(fit, tracks = NULL){

  if (!inherits(fit, "fitLangevin")) {
    stop("Input 'fit' must be a fitLangevin object.")
  }
  if (is.null(fit$osa)) {
    stop("The fit object does not contain OSA residuals. Make sure to run fitLangevin with calcOSA = TRUE.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting residuals. Please install it.")
  }

  if (is.null(tracks)) {
    return(generate_plots(fit$osa))
  } else {
    if (length(tracks) == 1 && tracks[1] == "all") {
      tracks_to_plot <- unique(fit$osa$id)
    } else {
      tracks_to_plot <- tracks
    }

    out_list <- list()
    for (trk in tracks_to_plot) {
      sub_data <- fit$osa[fit$osa$id == trk, ]
      if (nrow(sub_data) > 0) {
        out_list[[as.character(trk)]] <- generate_plots(sub_data, paste0(" - ID: ", trk))
      } else warning("No data found for track ID: ", trk, ". Skipping this track.")
    }
    return(out_list)
  }
}

make_acf_plot <- function(x, title) {
  x <- stats::na.omit(x)
  if(length(x) < 2) return(NULL)

  acf_obj <- stats::acf(x, plot = FALSE)
  acf_df <- data.frame("lag" = acf_obj$lag[, 1, 1], "acf" = acf_obj$acf[, 1, 1])
  clim <- stats::qnorm(0.975) / sqrt(acf_obj$n.used)

  ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = c(0, -clim, clim),
                        color = c("black", "blue", "blue"),
                        linetype = c("solid", "dashed", "dashed")) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::labs(title = title, x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
}

make_qq_norm <- function(x, title) {
  x <- sort(stats::na.omit(x))
  n <- length(x)
  if(n < 2) return(NULL)

  p <- stats::ppoints(n)
  z <- stats::qnorm(p)

  Qx <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
  Qz <- stats::qnorm(c(0.25, 0.75))
  slope <- diff(Qx) / diff(Qz)
  int <- Qx[1] - slope * Qz[1]

  se <- slope / stats::dnorm(z) * sqrt(p * (1 - p) / n)
  fit_line <- int + slope * z
  upper <- fit_line + 1.96 * se
  lower <- fit_line - 1.96 * se

  outlier <- x < lower | x > upper

  df <- data.frame("sample" = x, "theoretical" = z, "lower" = lower, "upper" = upper, "outlier" = outlier)

  ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
    ggplot2::geom_abline(intercept = int, slope = slope, color = "red", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = outlier), alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "red")) +
    ggplot2::labs(title = title, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

make_qq_chisq <- function(x, title) {
  x <- sort(stats::na.omit(x))
  n <- length(x)
  if(n < 2) return(NULL)

  p <- stats::ppoints(n)
  z <- stats::qchisq(p, df = 2)

  slope <- 1
  int <- 0

  se <- slope / stats::dchisq(z, df = 2) * sqrt(p * (1 - p) / n)
  fit_line <- int + slope * z
  upper <- fit_line + 1.96 * se
  lower <- fit_line - 1.96 * se

  outlier <- x < lower | x > upper

  df <- data.frame("sample" = x, "theoretical" = z, "lower" = lower, "upper" = upper, "outlier" = outlier)

  ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
    ggplot2::geom_abline(intercept = int, slope = slope, color = "red", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = outlier), alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "red")) +
    ggplot2::labs(title = title, x = "Theoretical Quantiles (Chi-sq, df=2)", y = "Sample Squared Mahalanobis Distance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

generate_plots <- function(data_subset, title_suffix = "") {
  res_x_full <- data_subset$residual.x[!is.na(data_subset$residual.x)]
  res_y_full <- data_subset$residual.y[!is.na(data_subset$residual.y)]

  valid_idx <- which(!is.na(data_subset$residual.x) & !is.na(data_subset$residual.y))
  res_x_pair <- data_subset$residual.x[valid_idx]
  res_y_pair <- data_subset$residual.y[valid_idx]
  D2 <- res_x_pair^2 + res_y_pair^2

  p_qq_x <- make_qq_norm(res_x_full, paste0("Normal Q-Q Plot: x residuals", title_suffix))
  p_qq_y <- make_qq_norm(res_y_full, paste0("Normal Q-Q Plot: y residuals", title_suffix))
  p_acf_x <- make_acf_plot(res_x_full, paste0("ACF: x residuals", title_suffix))
  p_acf_y <- make_acf_plot(res_y_full, paste0("ACF: y residuals", title_suffix))

  p_qq_mah <- make_qq_chisq(D2, paste0("Mahalanobis Q-Q Plot", title_suffix))
  p_acf_mah <- make_acf_plot(D2, paste0("ACF: Mahalanobis Distance", title_suffix))

  return(list(
    qq_x = p_qq_x, qq_y = p_qq_y,
    acf_x = p_acf_x, acf_y = p_acf_y,
    qq_mah = p_qq_mah, acf_mah = p_acf_mah
  ))
}
