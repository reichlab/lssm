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
#' @param p autoregressive order
#'
#' @return arma fit object
#'
#' @export
new_arma <- function(
    param_estimates
) {
  arma <- structure(
    list(
      param_estimates = param_estimates,
      p = 1
    ),
    class = "arma"
  )

  return(arma)
}


#' Fit a quantile baseline model to historical disease incidence
#'
#' @param y numeric vector of response values
#' @param verbose logical: if TRUE, print output from estimation in Stan
#' @param ... other arguments are ignored
#'
#' @return arma fit object
#'
#' @export
fit_arma <- function(
    y,
    verbose = FALSE,
    ...) {
  stan_path <- file.path(
    find.package("lssm"),
    "stan_models",
    "AR1.stan"
  )

  model <- rstan::stan_model(stan_path)

  stan_data <- list(
    n = length(y),
    p = 1,
    y = array(data = y, dim = c(length(y), 1))
  )

  estimate <- rstan::optimizing(
    object = model,
    data = stan_data,
    verbose = verbose)

  return(new_arma(
    param_estimates = estimate$par
  ))
}


#' Predict future disease incidence starting from the end of the training data.
#'
#' @param arma_fit an arma fit object
#' @param post_pred logical; if TRUE, samples from the posterior predictive are
#' returned for time points with observations in the training set
#' @param horizon number of time steps forward to predict
#' @param forecast_representation string specifying approach to representing
#' forecast distributions, see documentation of represent_forecasts
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
  post_pred,
  horizon,
  forecast_representation,
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1e5) {
  if(horizon > 1) {
    stop("forecast horizon > 1 not yet supported")
  }

  stan_data <- list(
    n = length(newdata),
    p = 1,
    y = array(data = newdata, dim = c(length(newdata), 1)),
    steps_ahead = horizon,
    m = 1,
    phi_1 = arma_fit$param_estimates["phi_1"],
    var_zeta =  arma_fit$param_estimates["var_zeta"],
    a1 = arma_fit$param_estimates["a1[1]"]
  )
  dim(stan_data$a1) <- 1

  stan_path <- file.path(
    find.package("lssm"),
    "stan_models",
    "predict.stan"
  )

  model <- rstan::stan_model(stan_path)

  stan_output <- rstan::sampling(
    model,
    data = stan_data,
    algorithm = "Fixed_param",
    chains = 1, iter = 1, warmup = 0
  )

  raw_prediction <- rstan::extract(stan_output)$results

  named_dist_forecast <- purrr::map_dfr(
    seq_len(horizon),
    function(h) {
      data.frame(
        h = h,
        family = "norm",
        mean = raw_prediction[1, 1, 1],
        sd = sqrt(raw_prediction[1, 1, 2]),
        stringsAsFactors = FALSE
      )
    }
  )

  forecast <- represent_forecasts(
    named_dist_forecast = named_dist_forecast,
    forecast_representation = forecast_representation,
    quantile_levels = quantile_levels,
    nsim = nsim
  )

  return(forecast)
}
