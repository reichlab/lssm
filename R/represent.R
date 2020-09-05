#' Translate representation of forecast distributions to a standardized form
#'
#' @param named_dist_forecast data frame specifying the forecast, with at
#' minimum a "family" column. Any columns occurring before the "family"
#' are assumed to be identifiers of the forecast, such as a forecast horizon or
#' location; these columns are retained in the output. Any columns after the
#' "family" are assumed to be parameters of the forecast distribution, and are
#' used as part of calls to the q or d function for the family.
#' @param forecast_representation string specifying approach to representing
#' forecast distributions, see documentation of represent_forecasts
#' @param quantile_levels numeric vector of quantile levels to use for
#' forecast_representation = "quantile"
#' @param nsim integer number of samples to use for
#' forecast_representation = "sample"
#'
#' @return a representation of the forecast; the structure depends on
#' `forecast_representation`:
#' * "named_dist": the `named_dist_forecast` argument is returned
#' * "sample": a matrix is returned with number of rows equal to the number of
#' rows of `named_dist_forecast` and number of columns equal to `nsim`.
#' Each row contains `nsims` samples from the predictive distribution in the
#' corresponding row of `named_dist_forecast`.
#' * "quantile": a matrix is returned with number of rows equal to the number
#' of rows of `named_dist_forecast` and number of columns equal to the length
#' of `quantile_levels`. Each row contains quantiles of the predictive
#' distribution in the corresponding row of `named_dist_forecast`.
#'
#' @export
represent_forecasts <- function(
  named_dist_forecast,
  forecast_representation = c("named_dist", "sample", "quantile"),
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1e5
) {
  forecast_representation <- match.arg(
    forecast_representation,
    choices = c("named_dist", "sample", "quantile"))

  if (forecast_representation == "named_dist") {
    forecast <- named_dist_forecast
  } else if (forecast_representation == "sample") {
    fun_name <- paste0("r", named_dist_forecast$family[1])
    family_col_ind <- which(colnames(named_dist_forecast) == "family")
    param_cols <- family_col_ind +
      seq_len(ncol(named_dist_forecast) - family_col_ind)
    forecast <- matrix(
      NA_real_,
      nrow = nrow(named_dist_forecast),
      ncol = nsim)
    colnames(forecast) <- paste0("sample_", seq_len(nsim))

    for (h in seq_len(nrow(named_dist_forecast))) {
      call_args <- as.list(named_dist_forecast[h, param_cols])
      call_args$n <- nsim
      forecast[h, ] <- do.call(what = fun_name, args = call_args)
    }
  } else if (forecast_representation == "quantile") {
    fun_name <- paste0("q", named_dist_forecast$family[1])
    family_col_ind <- which(colnames(named_dist_forecast) == "family")
    param_cols <- family_col_ind +
      seq_len(ncol(named_dist_forecast) - family_col_ind)
    forecast <- matrix(
      NA_real_,
      nrow = nrow(named_dist_forecast),
      ncol = length(quantile_levels))
    colnames(forecast) <- paste0("quantile_", quantile_levels)

    for (h in seq_len(nrow(named_dist_forecast))) {
      call_args <- as.list(named_dist_forecast[h, param_cols])
      call_args$p <- quantile_levels
      forecast[h, ] <- do.call(what = fun_name, args = call_args)
    }
  }

  return(forecast)
}
