## utility functions for lssm models -- mainly, transformations and differencing

#' Do an initial transformation of time series data.
#'
#' @param y a univariate time series or numeric vector.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", or "none".
#' @param transform_offset numeric offset used before the Box-Cox and log
#'   transformations; the offset is added to all observations before
#'   transforming.  Default value of 0.5 allows us to use the Box-Cox and log
#'   transforms (which require positive inputs) in case of observations of 0,
#'   and also ensures that the reverse transformed values will always be at
#'   least -0.5, so that they round up to non-negative values.
#' @param bc_lambda (required if transformation is "box-cox").  Parameter for
#'   the Box-Cox transformation.  bc_lambda should be the result of a
#'   call to car::powerTransform(y + transform_offset, family = "bcPower")
#'
#' @return a transformed object of the same class as y
#'
#' @export
do_initial_transform <- function(
  y,
  transformation,
  transform_offset = 0.5,
  bc_lambda) {
  if (identical(transformation, "log")) {
    transformed_y <- log(y + transform_offset)
  } else if (identical(transformation, "box-cox")) {
    if (missing(bc_lambda)) {
      stop("bc_lambda must be provided if transformation is 'box-cox'")
    }
    transformed_y <- car::bcPower(
      U = y + transform_offset,
      lambda = bc_lambda)
  } else if (identical(transformation, "none")) {
    transformed_y <- y
  } else {
    stop('Invalid transformation: must be one of "box-cox", "log", or "none".')
  }

  return(transformed_y)
}

#' Invert an initial transformation of time series data.
#'
#' @param y a univariate time series or numeric vector.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", or "none".
#' @param transform_offset numeric offset used before the Box-Cox and log
#'   transformations; the offset is added to all observations before
#'   transforming.  Default value of 0.5 allows us to use the Box-Cox and log
#'   transforms (which require positive inputs) in case of observations of 0,
#'   and also ensures that the reverse transformed values will always be at
#'   least -0.5, so that they round up to non-negative values.
#' @param bc_lambda (required if transformation is "box-cox").  Parameter for
#'   the Box-Cox transformation.  bc_lambda should be the result of a
#'   call to car::powerTransform(y + transform_offset, family = "bcPower")
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_initial_transform <- function(
  y,
  transformation,
  transform_offset = 0.5,
  bc_lambda) {
  if (identical(transformation, "log")) {
    detransformed_y <- exp(y) - transform_offset
  } else if (identical(transformation, "box-cox")) {
    if (missing(bc_lambda)) {
      stop("bc_lambda must be provided if transformation is 'box-cox'")
    }
    detransformed_y <- invert_bc_transform(b = y,
      lambda = bc_lambda,
      gamma = transform_offset)
  } else if (identical(transformation, "none")) {
    detransformed_y <- y
  } else {
    stop("Invalid transformation: must be one of 'box-cox', 'log', or 'none'.")
  }

  return(detransformed_y)
}

#' Invert a Box-Cox transformation.  See car::bcPower for the original
#' transformation.
#'
#' @param b a univariate numeric vector.
#' @param lambda exponent for Box-Cox transformation
#' @param gamma offset for Box-Cox transformation
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_bc_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) subtract offset gamma

  ## 1) undo box-cox
  if (abs(lambda) <= 1e-10) {
    z <- exp(b)
  } else {
    z <- (lambda * b + 1)^(1 / lambda)
  }

  ## 2) Subtract gamma
  return(z - gamma)
}

#' Invert an initial transformation of time series data.
#'
#' @param y a univariate time series or numeric vector.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", or "none".
#' @param transform_offset numeric offset used before the Box-Cox and log
#'   transformations; the offset is added to all observations before
#'   transforming.  Default value of 0.5 allows us to use the Box-Cox and log
#'   transforms (which require positive inputs) in case of observations of 0,
#'   and also ensures that the reverse transformed values will always be at
#'   least -0.5, so that they round up to non-negative values.
#' @param bc_lambda (required if transformation is "box-cox").  Parameter for
#'   the Box-Cox transformation.  bc_lambda should be the result of a
#'   call to car::powerTransform(y + transform_offset, family = "bcPower")
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_initial_transform_probabilistic <- function(
  y,
  transformation,
  transform_offset = 0.5,
  bc_lambda) {
  detransformed_y <- y
  if (identical(transformation, "log")) {
    detransformed_y$family <- "logmvnorm"
    detransformed_y$offset <- transform_offset
  } else if (identical(transformation, "box-cox")) {
    detransformed_y$family <- "bcmvnorm"
    detransformed_y$lambda <- bc_lambda
    detransformed_y$offset <- transform_offset
  } else if (!identical(transformation, "none")) {
    stop("Invalid transformation: must be one of 'box-cox', 'log', or 'none'.")
  }
  
  return(detransformed_y)
}


#' Do first-order and seasonal differencing (go from original time series
#' to differenced time series).
#'
#' @param y a univariate time series or numeric vector.
#' @param d order of first differencing
#' @param D order of seasonal differencing
#' @param frequency frequency of time series.  Must be provided if y is not
#'   of class "ts" and D > 0.  See the help for stats::ts for more.
#'
#' @return a differenced time series object (of class 'ts'),
#'   padded with leading NAs.
#'
#' @export
do_difference <- function(y, d = 0, D = 0, frequency = 1) {
  # first differencing
  for (i in seq_len(d)) {
#    y <- ts(
    y <- c(NA,
        y[seq(from = 1 + 1, to = length(y))] -
          y[seq(from = 1, to = length(y) - 1)])#,
#    frequency = frequency)
  }

  # seasonal differencing
  if (D > 0 && frequency < 2) {
    stop(paste0("It doesn't make sense to do seasonal differencing with a ",
                "time series frequency of 1."))
  }
  for (i in seq_len(D)) {
#    y <- ts(
    y <- c(rep(NA, frequency),
        y[seq(from = frequency + 1, to = length(y))] -
          y[seq(from = 1, to = length(y) - frequency)])#,
#      frequency = frequency)
  }

  return(y)
}


#' Invert first-order and seasonal differencing (go from seasonally differenced
#' time series to original time series).  In this function, all values of y and
#' dy are fixed numeric values.
#'
#' @param dy a first-order and/or seasonally differenced univariate time series
#'   with values like y_{t} - y_{t - ts_frequency}
#' @param y a univariate time series or numeric vector with values like
#'   y_{t - ts_frequency}.
#' @param d order of first differencing
#' @param D order of seasonal differencing
#' @param frequency frequency of time series.  Must be provided if y is not
#'   of class "ts" and D > 0.  See the help for stats::ts for more.
#'
#' @details y may have longer length than dy.  It is assumed that dy "starts"
#'   one time index after y "ends": that is, if y is of length T, d = 0, and
#'   D = 1 then dy[1] = y[T + 1] - y[T + 1 - ts_frequency]
#'
#' @return a time series object (of class 'ts')
#'
#' @export
invert_difference_deterministic <- function(dy, y, d, D, frequency) {
  for (i in seq_len(d)) {
    y_dm1 <- do_difference(y, d = d - i, D = D, frequency = frequency)
    dy_full <- c(y_dm1, dy)
    for (t in seq_len(length(dy))) {
      dy_full[length(y_dm1) + t] <- dy_full[length(y_dm1) + t - 1] +
        dy_full[length(y_dm1) + t]
    }
    dy <- dy_full[length(y_dm1) + seq_along(dy)]
  }

  for (i in seq_len(D)) {
    y_dm1 <- do_difference(y, d = 0, D = D - i, frequency = frequency)
    dy_full <- c(y_dm1, dy)
    for (t in seq_len(length(dy))) {
      dy_full[length(y_dm1) + t] <- dy_full[length(y_dm1) + t - frequency] +
        dy_full[length(y_dm1) + t]
    }
    dy <- dy_full[length(y_dm1) + seq_along(dy)]
  }

  return(dy)
}


#' Invert first-order and seasonal differencing (go from seasonally differenced
#' time series to original time series).  In this function, all values of y are
#' fixed numeric values and dy must be a representation of the joint
#' distribution over all relevant time points with a multivariate normal
#' distribution.  This should be in the form of a tibble with columns mean
#' containing the mean and sigma containing the covariance matrix.
#'
#' @param dy a first-order and/or seasonally differenced univariate time series
#'   with values like y_{t} - y_{t - ts_frequency}
#' @param y a univariate time series or numeric vector with values like
#'   y_{t - ts_frequency}.
#' @param d order of first differencing
#' @param D order of seasonal differencing
#' @param frequency frequency of time series.  Must be provided if y is not
#'   of class "ts" and D > 0.  See the help for stats::ts for more.
#'
#' @details y may have longer length than dy.  It is assumed that dy "starts"
#'   one time index after y "ends": that is, if y is of length T, d = 0, and
#'   D = 1 then dy[1] = y[T + 1] - y[T + 1 - ts_frequency]
#'
#' @return a tibble with all the same entries as the argument dy, but with
#' the mean and covariance of the original undifferenced data at those time
#' points
#'
#' @export
invert_difference_probabilistic <- function(dy, y, d, D, frequency) {
  if(d > 0) {
    # set up matrix encoding un-differencing operation
    A <- matrix(0, nrow = length(dy$mean[[1]]), ncol = length(dy$mean[[1]]) + 1)
    for(i in seq_len(nrow(A))) {
      A[i, seq_len(i + 1)] <- 1L
    }
    
    # augmented covariance with 0s for observed quantities
    sigma_full <- matrix(
      0,
      nrow = length(dy$mean[[1]]) + 1,
      ncol = length(dy$mean[[1]]) + 1
    )
    
    for (i in seq_len(d)) {
      y_dm1 <- do_difference(y, d = d - i, D = D, frequency = frequency)
      mean_full <- c(tail(y_dm1, 1), dy$mean[[1]])
      sigma_full[1 + seq_along(dy$mean[[1]]), 1 + seq_along(dy$mean[[1]])] <-
        dy$sigma[[1]]
      dy$mean[[1]] <- A %*% mean_full
      dy$sigma[[1]] <- emulator::quad.tform(sigma_full, A)
    }
  }
  
  if(D > 0) {
    # set up matrix encoding un-differencing operation
    A <- matrix(0, nrow = length(dy$mean[[1]]), ncol = length(dy$mean[[1]]) + frequency)
    for(i in seq_len(nrow(A))) {
      A[i, seq(from = i %% frequency, by = frequency, to = i + frequency)] <- 1L
    }
    
    # augmented covariance with 0s for observed quantities
    sigma_full <- matrix(
      0,
      nrow = length(dy$mean[[1]]) + frequency,
      ncol = length(dy$mean[[1]]) + frequency
    )
    
    for (i in seq_len(D)) {
      y_dm1 <- do_difference(y, d = 0, D = D - i, frequency = frequency)
      mean_full <- c(tail(y_dm1, frequency), dy$mean[[1]])
      sigma_full[frequency + seq_along(dy$mean[[1]]),
                 frequency + seq_along(dy$mean[[1]])] <- dy$sigma[[1]]
      dy$mean[[1]] <- A %*% mean_full
      dy$sigma[[1]] <- emulator::quad.tform(sigma_full, A)
    }
  }
  
  return(dy)
}



#' Remove leading values that are infinite or missing, and replace all internal
#' sequences of infinite or missing values by linearly interpolating between
#' the surrounding observed (finite) values.
#'
#' @param y a univariate time series or numeric vector
#'
#' @return a vector of the same class as y
#'
#' @export
interpolate_and_clean_missing <- function(y) {
  if (any(is.na(y) | is.infinite(y)) && !is.na(tail(y, 1))) {
    ## drop leading NAs
    y <- drop_leading_nas(y)

    ## interpolate internal NAs
    while (any(is.na(y) | is.infinite(y))) {
      na_rle <- rle(is.na(y) | is.infinite(y))
      na_run_start_ind <- na_rle$lengths[1] + 1
      na_run_end_ind <- na_run_start_ind + na_rle$lengths[2] - 1
      y[na_run_start_ind:na_run_end_ind] <-
        approx(
          x = c(na_run_start_ind - 1, na_run_end_ind + 1),
          y = y[c(na_run_start_ind - 1, na_run_end_ind + 1)],
          xout = na_run_start_ind:na_run_end_ind,
          method = "linear"
          )$y
    }
  }

  return(y)
}


#' Remove leading values that are infinite or missing
#'
#' @param y a univariate time series or numeric vector
#'
#' @return a vector of the same class as y
#'
#' @export
drop_leading_nas <- function(y) {
  ## drop leading NAs
  if (is.na(y[1]) | is.infinite(y[1])) {
    num_leading_nas <- rle(is.na(y) | is.infinite(y))$lengths[1]
    y <- y[- seq_len(num_leading_nas)]
  }

  return(y)
}
