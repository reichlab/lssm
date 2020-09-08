#' Create a sequence of integers
#'
#' @param from starting value for sequence
#' @param to ending value for sequence
#' @param max_len maximum length of sequence
#' 
#' @return sequence of integers starting at from, ending at to, and of length
#' at most max_len, with a preference for smaller integers
int_seq <- function(from, to, max_len) {
  len <- min(max_len, to - from + 1)
  seq(from = from, to = to, len = len) %>%
    floor()
}

#' Define grid of tuning parameters for the arma model
#' 
#' @param x 
#' @param y 
#' @param len 
#' @param search 
#' @param transform 
#' @param transform_offset 
#' @param max_d 
#' @param max_D 
#' @param max_p_ar 
#' @param max_q_ma 
#' @param min_order
#' @param max_order 
#' 
#' @return data frame with columns model, transform, transform_offset, d, D,
#' p_ar, q_ma
#' 
#' @export
arma_param_grid <- function(
  x,
  y,
  len = NULL,
  search = "grid",
  transformation = "box-cox",
  transform_offset = if (any(x <= 0)) { -1 * min(x) + 0.49} else { 0.0 },
  max_d = 2,
  max_D = 1,
  max_p_ar = 5,
  max_q_ma = 5,
  min_order = 1,
  max_order = 5
) {
  library(lssm)
  
  # Only do seasonal differencing if x has a seasonal period
  if(frequency(x) < 2) {
    max_D <- 0
  }
  
  # For grid search, all combinations of d, D, p_ar, and q_ma satisfying
  # p_ar + q_ma <= max_order
  out <- expand.grid(
    model = "arma",
    transformation = transformation,
    transform_offset = transform_offset,
    d = int_seq(from = 0, to = max_d, max_len = len),
    D = int_seq(from = 0, to = max_D, max_len = len),
    p_ar = int_seq(from = 0, to = max_p_ar, max_len = len),
    q_ma = int_seq(from = 0, to = max_q_ma, max_len = len),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      p_ar + q_ma >= min_order,
      p_ar + q_ma <= max_order
    )
  
  # If data contains non-positive values, remove grid combinations where the
  # transformation is log or box-cox and the transform_offset is less than or
  # equal to the smallest data value
  min_data <- min(c(x, y))
  if(min_data <= 0) {
    out <- out %>%
      dplyr::filter(
        !(transformation %in% c("log", "box-cox") &
            transform_offset <= min_data)
      )
  }
  
  if (search == "random") {
    # select a random set of rows from the grid defined above
    out <- out %>%
      dplyr::sample_n(size = len)
  }
  
  return(out)
}


#' Wrapper around fit_lssm for use with caret::train
#' 
#' @param x time series data to fit to
#' @param y ignored
#' @param param dataframe of one row of arguments to fit_lssm
#' @param ... other arguments are ignored
#' 
#' @return numeric vector of predictive medians with attributes:
#'  * family is a string with the parametric family, e.g. "norm"
#'  * other attributes are names of parameters for the parametric family
fit_lssm_caret_wrapper <- function(
  x,
  y,
  param,
  ts_frequency = 1,
  verbose = FALSE,
  ...
) {
  param <- as.list(param)
  param$y <- x[, 1]
  param$ts_frequency <- ts_frequency
  param$verbose <- verbose
  
  do.call(fit_lssm, param)
}


#' Wrapper around predict.lssm for use with caret::train
#' 
#' @param model_fit model fit object of class lssm
#' @param newdata data from which to generate predictions
#' @param ... other arguments are ignored
#' 
#' @return numeric vector of predictive medians with attributes:
#'  * family is a string with the parametric family, e.g. "norm"
#'  * other attributes are names of parameters for the parametric family
predict_lssm_caret_wrapper <- function(
  modelFit,
  newdata,
  ...
) {
  named_dist_forecast <- predict(
    modelFit,
    horizon = nrow(newdata),
    forecast_representation = "named_dist")
  
  median_forecast <- predict(
    modelFit,
    horizon = nrow(newdata),
    forecast_representation = "quantile",
    quantile_levels = 0.5
  )[, 1]
  
  for(cname in colnames(named_dist_forecast)) {
    attr(median_forecast, cname) <- named_dist_forecast[[cname]]
  }
  
  return(median_forecast)
}

#' Calculate the log score for forecasts in the format required by caret
#' 
#' @param data data frame with columns obs and pred.  pred must be in the form
#' returned by predict_lssm_caret_wrapper
#' @param ... other arguments are ignored
#' 
#' @export
log_score_summary <- function(
  data,
  ...) {
  pred_attrs <- attributes(data$pred)
  
  if(is.null(pred_attrs)) {
    # this condition will run if log_score_summary is called with data$pred not
    # produced by predict_lssm_caret_wrapper.  This happens once inside of
    # caret::train when it checks that the summary function behaves correctly
    # so we have to return something rather than throw an error.
    return(c("log_score" = 1.0))
  }
  
  dfun <- paste0("d", pred_attrs$family[1])
  call_args <- pred_attrs[!(names(pred_attrs) %in% c("family", "h"))]
  call_args$x <- data$obs
  call_args$log <- TRUE
  
  return(c(
    "log_score" = mean(do.call(dfun, call_args))
  ))
}


#' @export
lssm_arma_caret <- list(
  library = "lssm",
  type = "Regression",
  parameters = dplyr::bind_rows(
    data.frame(
      parameter = "model",
      class = "character",
      label = "model",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "transformation",
      class = "character",
      label = "transformation",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "transform_offset",
      class = "numeric",
      label = "transform_offset",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "d",
      class = "numeric",
      label = "d",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "D",
      class = "numeric",
      label = "D",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "p_ar",
      class = "numeric",
      label = "p_ar",
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = "q_ma",
      class = "numeric",
      label = "q_ma",
      stringsAsFactors = FALSE
    )
  ),
  grid = arma_param_grid,
  fit = fit_lssm_caret_wrapper,
  predict = predict_lssm_caret_wrapper,
  prob = NULL#,
#  sort = sort_arma_params
)
