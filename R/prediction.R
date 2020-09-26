## functions for prediction from linear state space models

#' Generate predictions from a linear state space model fit
#'
#' This function handles any transformations and differencing that were done
#' before calling the model-specific predict method.
#'
#' @param object a model fit of class "lssm", as returned by fit_lssm
#' @param newdata numeric vector of new data to simulate forward from
#' @param forecast_representation string specifying approach to representing
#' forecast distributions. One of "named_dist", "sample", or "quantile";
#' see documentation of `represent_forecasts` for more detail.
#' @param joint logical; if TRUE, named distribution representation tracks the
#' joint distribution of y_t over all horizons.  if FALSE, we track only the
#' marginal distributions at each horizon.
#' @param quantile_levels numeric vector of quantile levels to use for
#' forecast_representation = "quantile"
#' @param nsim integer number of samples to use for
#' forecast_representation = "sample"
#' @param seed either `NULL` or an integer that will be used in a call to
#' `set.seed` before simulating the response vectors.  If set, the value is
#' saved as the "seed" attribute of the returned value.  The default, `NULL`,
#' will not change the random generator state, and return `.Random.seed`
#' as the "seed" attribute
#' @param horizon number of time steps forwards to simulate
#' @param ... other arguments passed on to model-specific predict methods
#'
#' @return an nsim by horizon matrix with simulated values
#'
#' @export
predict.lssm <- function(
  lssm_fit,
  newdata,
  forecast_representation,
  joint = TRUE,
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1,
  seed = NULL,
  horizon = 1,
  ...
) {
  forecast_representation <- match.arg(
    forecast_representation,
    choices = c("named_dist", "sample", "quantile"))
  
  if (is.null(seed)) {
    seed <- .Random.seed
  } else {
    set.seed(seed)
  }

  if (missing(newdata)) {
    newdata <- attr(lssm_fit, "lssm_y")
  }

  transformation <- attr(lssm_fit, "lssm_transformation")
  transform_offset <- attr(lssm_fit, "lssm_transform_offset")
  d <- attr(lssm_fit, "lssm_d")
  D <- attr(lssm_fit, "lssm_D")
  ts_frequency <- attr(lssm_fit, "lssm_ts_frequency")

  # Initial transformation, if necessary
  if (identical(transformation, "box-cox")) {
    bc_lambda <- attr(lssm_fit, "lssm_bc_lambda")
  } else {
    bc_lambda <- NULL
  }
  transformed_y <- do_initial_transform(
    y = newdata,
    transformation = transformation,
    transform_offset = transform_offset,
    bc_lambda = bc_lambda)

  # Initial differencing, if necessary
  differenced_y <- do_difference(transformed_y, d = d, D = D,
    frequency = ts_frequency)

  # Forecasts on the transformed and differenced scale
  raw_forecast <- NextMethod(newdata = differenced_y)

  # Forecasts are of seasonally differenced transformed time series
  # Get to forecasts for originally observed time series ("orig") by
  # inverting the differencing and transformation operations
  orig_forecast <- raw_forecast
  if(forecast_representation %in% c("sample", "quantile")) {
    if(any(grepl("quantile", colnames(orig_forecast))) &&
       nrow(orig_forecast) > 1) {
      stop(paste0("Differencing inversion for quantile forecasts at ",
        "horizons greater than 1 is not currently supported."))
    }
    
    for(i in seq_len(ncol(orig_forecast))) {
      orig_forecast[, i] <-
        invert_difference_deterministic(
          dy = orig_forecast[, i],
          y = transformed_y,
          d = d,
          D = D,
          frequency = ts_frequency)
  
      orig_forecast[, i] <-
        invert_initial_transform(
          y = orig_forecast[, i],
          transformation = transformation,
          transform_offset = transform_offset,
          bc_lambda = bc_lambda)
    }
  
    attr(orig_forecast, "seed") <- seed
  } else {
    if(nrow(orig_forecast) == 1 && orig_forecast$family == "mvnorm") {
      orig_forecast <- invert_difference_probabilistic(
        dy = orig_forecast,
        y = transformed_y,
        d = d,
        D = D,
        frequency = ts_frequency)
      
      orig_forecast <- invert_initial_transform_probabilistic(
        y = orig_forecast,
        transformation = transformation,
        transform_offset = transform_offset,
        bc_lambda = bc_lambda
      )
    } else {
      stop(paste0("Can only invert differencing with a joint forecast ",
        "distribution."))
    }
  }

  return(orig_forecast)
}
