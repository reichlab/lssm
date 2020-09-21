#' Check if object is of class arma
#'
#' @param object an object that may be a arma object
#'
#' @return boolean; whether object is inherits arma class
#'
#' @export
is.arma <- function(object) {
  return(inherits(object, "arma"))
}


#' Create a arma object
#'
#' @param param_estimates model parameter estimates
#' @param p_ar integer order for autoregressive terms
#' @param q_ma integer order of moving average terms
#'
#' @return arma fit object
#'
#' @export
new_arma <- function(
    param_estimates,
    p_ar,
    q_ma
) {
  arma <- structure(
    list(
      param_estimates = param_estimates,
      p_ar = p_ar,
      q_ma = q_ma
    ),
    class = "arma"
  )

  return(arma)
}


#' Fit a quantile baseline model to historical disease incidence
#'
#' @param y numeric vector of response values
#' @param p_ar integer order for autoregressive terms
#' @param q_ma integer order of moving average terms
#' @param verbose logical: if TRUE, print output from estimation in Stan
#' @param ... other arguments are ignored
#'
#' @return arma fit object
#'
#' @export
fit_arma <- function(
    y,
    p_ar,
    q_ma,
    verbose = FALSE,
    ...) {
  # drop leading NA's that may have resulted from differencing
  y <- drop_leading_nas(y)
  
  # input data for estimation
  stan_data <- list(
    n = length(y),
    p = 1,
    y = array(data = y, dim = c(length(y), 1)),
    p_ar = p_ar,
    q_ma = q_ma
  )

  # initial values for estimation
  # if (p_ar > 0) {
  #   # initialize to parameters for an AR(1) model
  #   acf_1 <- as.numeric(acf(y, lag.max = 1, plot = FALSE)["1"]$acf)
  # 
  #   init_par <- list(
  #     phi = c(acf_1, rep(0.0, p_ar - 1)),
  #     theta = rep(0.0, q_ma),
  #     var_zeta = var(y) * (1 - acf_1^2),
  #     a1 = c(y[1], rep(0.0, max(p_ar, q_ma + 1) - 1))
  #   )
  #   dim(init_par$a1) = length(init_par$a1)
  #   dim(init_par$phi) = length(init_par$phi)
  #   dim(init_par$theta) = length(init_par$theta)
  # } else if(q_ma > 0) {
  #   # initialize to something near the parameters for an MA(1) model
  #   acf_1 <- as.numeric(acf(y, lag.max = 1, plot = FALSE)["1"]$acf)
  #   if(1 - 4 * acf_1^2 < 0) {
  #     init_theta <- 1 / 2 * acf_1
  #   } else {
  #     init_theta <- (1 - sqrt(1 - 4 * acf_1^2)) / 2 * acf_1
  #   }
  # 
  #   init_par <- list(
  #     phi = rep(0.0, p_ar),
  #     theta = c(
  #       init_theta,
  #       rep(0.0, q_ma - 1)
  #     ),
  #     var_zeta = var(y) / (1 + init_theta^2),
  #     a1 = c(y[1], rep(0.0, max(p_ar, q_ma + 1) - 1))
  #   )
  #   dim(init_par$a1) = length(init_par$a1)
  #   dim(init_par$phi) = length(init_par$phi)
  #   dim(init_par$theta) = length(init_par$theta)
  # 
  # }# else {
    # initialize variance
    init_par <- list(
      phi = rep(0.0, p_ar),
      theta = rep(0.0, q_ma),
      var_zeta = var(y),
      a1 = c(y[1], rep(0.0, max(p_ar, q_ma + 1) - 1))#y[1]
    )
    dim(init_par$a1) = length(init_par$a1)
    dim(init_par$phi) = length(init_par$phi)
    dim(init_par$theta) = length(init_par$theta)
#  }

  estimate <- rstan::optimizing(
    object = stanmodels$ARMA_model,
    data = stan_data,
    init = init_par,
    verbose = verbose)

  return(new_arma(
    param_estimates = estimate$par,
    p_ar = p_ar,
    q_ma = q_ma
  ))
}


#' Predict future disease incidence starting from the end of the training data.
#'
#' @param arma_fit an arma fit object
#' @param horizon number of time steps forward to predict
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
#' @param ... other arguments are ignored
#'
#' @return representation of forecasts
#'
#' @export
predict.arma <- function(
  arma_fit,
  newdata,
  horizon,
  forecast_representation,
  joint = TRUE,
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1e5) {
  # drop leading NA's that may have resulted from differencing
  newdata <- drop_leading_nas(newdata)
  
  if (horizon > 1) {
    if(!joint) {
      stop("forecast horizon > 1 not yet supported")
    }
  }

  stan_data <- list(
    n = length(newdata),
    p = 1,
    y = array(data = newdata, dim = c(length(newdata), 1)),
    horizon = horizon,
    joint = as.integer(joint),
    p_ar = arma_fit$p_ar,
    q_ma = arma_fit$q_ma,
    r = max(arma_fit$p_ar, arma_fit$q_ma + 1),
    m = max(arma_fit$p_ar, arma_fit$q_ma + 1),
    phi = arma_fit$param_estimates[
      grepl("^phi", names(arma_fit$param_estimates))],
    theta = arma_fit$param_estimates[
      grepl("^theta", names(arma_fit$param_estimates))],
    var_zeta =  arma_fit$param_estimates["var_zeta"],
    a1 = arma_fit$param_estimates[
      grepl("^a1", names(arma_fit$param_estimates))]
  )
  dim(stan_data$phi) <- length(stan_data$phi)
  dim(stan_data$theta) <- length(stan_data$theta)
  dim(stan_data$a1) <- length(stan_data$a1)

  stan_output <- rstan::sampling(
    stanmodels$ARMA_predict,
    data = stan_data,
    algorithm = "Fixed_param",
    chains = 1, iter = 1, warmup = 0
  )

  raw_prediction <- rstan::extract(stan_output)$forecasts

  if(joint) {
    named_dist_forecast <- tibble(
      h = list(seq_len(horizon)),
      family = "mvnorm",
      mean = list(raw_prediction[1, 1, , 1]),
      sigma = list(
        raw_prediction[1, 1, , seq(from = 2, to = dim(raw_prediction)[4])]
      )
    )
  } else {
    named_dist_forecast <- data.frame(
      h = seq_len(horizon),
      family = "norm",
      mean = raw_prediction[1, 1, 1],
      sd = sqrt(raw_prediction[1, 1, 2]),
      stringsAsFactors = FALSE
    )
  }

  forecast <- represent_forecasts(
    named_dist_forecast = named_dist_forecast,
    forecast_representation = forecast_representation,
    quantile_levels = quantile_levels,
    nsim = nsim
  )

  return(forecast)
}
