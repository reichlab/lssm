#' Check if object is of class sarima_lssm
#'
#' @param object an object that may be a sarima_lssm object
#'
#' @return boolean; whether object is inherits sarima_lssm class
#'
#' @export
is.sarima_lssm <- function(object) {
  return(inherits(object, "sarima_lssm"))
}


#' Create a sarima_lssm object
#'
#' @param param_estimates model parameter estimates
#' @param include_intercept logical; if TRUE, include intercept
#' @param p_ar integer order for non-seasonal autoregressive terms
#' @param q_ma integer order of non-seasonal moving average terms
#' @param P_ar integer order for seasonal autoregressive terms
#' @param Q_ma integer order of seasonal moving average terms
#' @param ts_frequency integer frequency of time series
#'
#' @return sarima_lssm fit object
#'
#' @export
new_sarima_lssm <- function(
    param_estimates,
    include_intercept,
    p_ar,
    q_ma,
    P_ar,
    Q_ma,
    ts_frequency
) {
  sarima_lssm <- structure(
    list(
      param_estimates = param_estimates,
      include_intercept = include_intercept,
      p_ar = p_ar,
      q_ma = q_ma,
      P_ar = P_ar,
      Q_ma = Q_ma,
      ts_frequency = ts_frequency
    ),
    class = "sarima_lssm"
  )

  return(sarima_lssm)
}


#' Fit a quantile baseline model to historical disease incidence
#'
#' @param y numeric vector of response values
#' @param include_intercept logical: if TRUE, include an intercept
#' @param p_ar integer order for non-seasonal autoregressive terms
#' @param q_ma integer order of non-seasonal moving average terms
#' @param P_ar integer order for seasonal autoregressive terms
#' @param Q_ma integer order of seasonal moving average terms
#' @param ts_frequency integer frequency of time series
#' @param init_par initial values for model parameters in estimation
#' @param verbose logical: if TRUE, print output from estimation in Stan
#' @param ... other arguments are ignored
#'
#' @return sarima_lssm fit object
#'
#' @export
fit_sarima_lssm <- function(
    y,
    include_intercept,
    p_ar,
    q_ma,
    P_ar,
    Q_ma,
    ts_frequency,
    init_par,
    verbose = FALSE,
    ...) {
  
  # drop leading NA's that may have resulted from differencing
  y <- drop_leading_nas(y)
  
  # input data for estimation
  stan_data <- list(
    n = length(y),
    p = 1,
    y = array(data = y, dim = c(length(y), 1)),
    include_intercept = as.integer(include_intercept),
    p_ar = p_ar,
    q_ma = q_ma,
    P_ar = P_ar,
    Q_ma = Q_ma,
    ts_frequency = ts_frequency
  )

  # initial values for parameter estimation
  if(missing(init_par) || is.null(init_par)) {
    # initial parameter values not provided. We use the following settings:
    # - 0 for all model coefficients
    # - var(y) for the variance of the MA or error terms
    # - [y[1], 0, ..., 0] for the initial state
    init_par <- list(
      phi_0 = rep(0.0, include_intercept),
      phi = rep(0.0, p_ar),
      theta = rep(0.0, q_ma),
      phi_seasonal = rep(0.0, P_ar),
      theta_seasonal = rep(0.0, Q_ma),
      var_zeta = var(y),
      a1 = c(y[1], rep(0.0, max(p_ar + P_ar * ts_frequency,
                                q_ma + Q_ma * ts_frequency + 1) - 1))
    )
  } else {
    if(is.list(init_par)) {
      # validate correct parameters were provided -- names and dimensions
      stop("initializing parameters with a list is not yet supported.")
    } else {
      init_par <- list(
        phi_0 = init_pars[grepl("^phi_0\\[", names(init_par))],
        phi = init_par[grepl("^phi\\[", names(init_par))],
        phi_seasonal = init_par[grepl("^phi_seasonal\\[", names(init_par))],
        theta = init_par[grepl("^theta\\[", names(init_par))],
        theta_seasonal = init_par[grepl("^theta_seasonal\\[", names(init_par))],
        var_zeta =  init_par["var_zeta"],
        a1 = init_par[grepl("^a1", names(init_par))]
      )
    }
  }
  dim(init_par$a1) <- length(init_par$a1)
  dim(init_par$phi_0) <- length(init_par$phi_0)
  dim(init_par$phi) <- length(init_par$phi)
  dim(init_par$theta) <- length(init_par$theta)
  dim(init_par$phi_seasonal) <- length(init_par$phi_seasonal)
  dim(init_par$theta_seasonal) <- length(init_par$theta_seasonal)
  
  estimate <- rstan::optimizing(
    object = stanmodels$SARIMA_model,
    data = stan_data,
    init = init_par,
    verbose = verbose)
  
  return(new_sarima_lssm(
    param_estimates = estimate$par,
    include_intercept = include_intercept,
    p_ar = p_ar,
    q_ma = q_ma,
    P_ar = P_ar,
    Q_ma = Q_ma,
    ts_frequency = ts_frequency
  ))
}


#' Predict future disease incidence starting from the end of the training data.
#'
#' @param sarima_lssm_fit an sarima_lssm fit object
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
predict.sarima_lssm <- function(
  sarima_lssm_fit,
  newdata,
  horizon,
  forecast_representation,
  joint = TRUE,
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1e5) {
  # extract time series frequency
  ts_frequency <- attr(sarima_lssm_fit, "lssm_ts_frequency")
  
  # drop leading NA's that may have resulted from differencing
  newdata <- drop_leading_nas(newdata)
  
  if (horizon > 1) {
    if(!joint) {
      stop("forecast horizon > 1 not yet supported")
    }
  }

  r <- max(
    sarima_lssm_fit$p_ar + sarima_lssm_fit$P_ar * ts_frequency,
    sarima_lssm_fit$q_ma + sarima_lssm_fit$Q_ma * ts_frequency + 1)
  stan_data <- list(
    n = length(newdata),
    p = 1,
    y = array(data = newdata, dim = c(length(newdata), 1)),
    horizon = horizon,
    joint = as.integer(joint),
    include_intercept = sarima_lssm_fit$include_intercept,
    p_ar = sarima_lssm_fit$p_ar,
    q_ma = sarima_lssm_fit$q_ma,
    P_ar = sarima_lssm_fit$P_ar,
    Q_ma = sarima_lssm_fit$Q_ma,
    ts_frequency = ts_frequency,
    r = r,
    m = r,
    phi_0 = sarima_lssm_fit$param_estimates[
      grepl("^phi_0\\[", names(sarima_lssm_fit$param_estimates))],
    phi = sarima_lssm_fit$param_estimates[
      grepl("^phi\\[", names(sarima_lssm_fit$param_estimates))],
    phi_seasonal = sarima_lssm_fit$param_estimates[
      grepl("^phi_seasonal\\[", names(sarima_lssm_fit$param_estimates))],
    theta = sarima_lssm_fit$param_estimates[
      grepl("^theta\\[", names(sarima_lssm_fit$param_estimates))],
    theta_seasonal = sarima_lssm_fit$param_estimates[
      grepl("^theta_seasonal\\[", names(sarima_lssm_fit$param_estimates))],
    var_zeta =  sarima_lssm_fit$param_estimates["var_zeta"],
    a1 = sarima_lssm_fit$param_estimates[
      grepl("^a1", names(sarima_lssm_fit$param_estimates))]
  )
  dim(stan_data$phi_0) <- length(stan_data$phi_0)
  dim(stan_data$phi) <- length(stan_data$phi)
  dim(stan_data$phi_seasonal) <- length(stan_data$phi_seasonal)
  dim(stan_data$theta) <- length(stan_data$theta)
  dim(stan_data$theta_seasonal) <- length(stan_data$theta_seasonal)
  dim(stan_data$a1) <- length(stan_data$a1)

  stan_output <- rstan::optimizing(
    stanmodels$SARIMA_predict,
    data = stan_data
  )
  
  raw_prediction <- stan_output$theta_tilde[,
    grepl("^forecasts\\[", colnames(stan_output$theta_tilde))]
  dim(raw_prediction) <- c(horizon, horizon + 1)
  
  if(joint) {
    named_dist_forecast <- tibble(
      h = list(seq_len(horizon)),
      family = "mvnorm",
      mean = list(raw_prediction[ , 1]),
      sigma = list(
        raw_prediction[, seq(from = 2, to = ncol(raw_prediction)), drop = FALSE]
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
