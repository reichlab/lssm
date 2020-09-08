## functions for prediction from linear state space models

#' Generate predictions from a linear state space model fit
#'
#' This function handles any transformations and differencing that were done
#' before calling the model-specific predict method.
#'
#' @param object a model fit of class "lssm", as returned by fit_lssm
#' @param newdata numeric vector of new data to simulate forward from
#' @param nsim number of sample trajectories to simulate
#' @param seed either `NULL` or an integer that will be used in a call to
#'   `set.seed` before simulating the response vectors.  If set, the value is
#'   saved as the "seed" attribute of the returned value.  The default, `NULL`,
#'   will not change the random generator state, and return `.Random.seed`
#'   as the "seed" attribute
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
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1,
  seed = NULL,
  horizon = 1,
  ...
) {
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
  if(is.matrix(orig_forecast)) {
    if(any(grepl("quantile", colnames(orig_forecast))) &&
       nrow(orig_forecast) > 1) {
      # Currently, just fail; we could do better.
      stop(paste0("Differencing inversion for quantile forecasts at ",
        "horizons grater than 1 is not currently supported."))
    }
    
    for(i in seq_len(ncol(orig_forecast))) {
      orig_forecast[, i] <-
        invert_difference(
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
  } else if (is.data.frame(orig_forecast)) {
    if(nrow(orig_forecast) > 1) {
      # Currently, just fail; we could do better. Doing so would need to
      # account for correlations across horizons if adding forecasts at
      # multiple horizons, to adjust variances
      stop(paste0("Differencing inversion for quantile forecasts at ",
                  "horizons grater than 1 is not currently supported."))
    }
    
    if(all(orig_forecast$family == "norm")) {
      orig_forecast$mean <- invert_difference(
        dy = orig_forecast$mean,
        y = transformed_y,
        d = d,
        D = D,
        frequency = ts_frequency)
    }
  }

  return(orig_forecast)
}
