## functions for prediction from linear state space models

#' Generate predictions from a linear state space model fit
#'
#' This function handles any transformations and differencing that were done
#' before calling the model-specific predict method.
#'
#' @param object a model fit of class "lssm", as returned by fit_lssm
#' @param nsim number of sample trajectories to simulate
#' @param seed either `NULL` or an integer that will be used in a call to
#'   `set.seed` before simulating the response vectors.  If set, the value is
#'   saved as the "seed" attribute of the returned value.  The default, `NULL`,
#'   will not change the random generator state, and return `.Random.seed`
#'   as the "seed" attribute
#' @param newdata numeric vector of new data to simulate forward from
#' @param horizon number of time steps forwards to simulate
#' @param ... other arguments passed on to model-specific predict methods
#'
#' @return an nsim by horizon matrix with simulated values
#'
#' @export
predict.lssm <- function(
  lssm_fit,
  nsim = 1,
  seed = NULL,
  newdata,
  horizon = 1,
  ...
) {
  if(is.null(seed)) {
    seed <- .Random.seed
  } else {
    set.seed(seed)
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

  # Sample trajectories on the transformed and differenced scale
  raw_trajectory_samples <- NextMethod(newdata = differenced_y)

  # Sampled trajectories are of seasonally differenced transformed time series
  # Get to trajectories for originally observed time series ("orig") by
  # adding inverting the differencing and transformation operations
  orig_trajectory_samples <- raw_trajectory_samples
  for(i in seq_len(nsim)) {
    orig_trajectory_samples[i, ] <-
      invert_difference(
        dy = raw_trajectory_samples[i, ],
        y = transformed_y,
        d = d,
        D = D,
        frequency = ts_frequency)

    orig_trajectory_samples[i, ] <-
      invert_initial_transform(
        y = orig_trajectory_samples[i, ],
        transformation = transformation,
        transform_offset = transform_offset,
        bc_lambda = bc_lambda)
  }

  attr(orig_trajectory_samples, "seed") <- seed

  return(orig_trajectory_samples)
}
