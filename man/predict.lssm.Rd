% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{predict.lssm}
\alias{predict.lssm}
\title{Generate predictions from a linear state space model fit}
\usage{
\method{predict}{lssm}(
  lssm_fit,
  newdata,
  forecast_representation,
  joint = TRUE,
  quantile_levels = c(0.025, 0.25, 0.5, 0.75, 0.975),
  nsim = 1,
  seed = NULL,
  horizon = 1,
  ...
)
}
\arguments{
\item{newdata}{numeric vector of new data to simulate forward from}

\item{forecast_representation}{string specifying approach to representing
forecast distributions. One of "named_dist", "sample", or "quantile";
see documentation of \code{represent_forecasts} for more detail.}

\item{joint}{logical; if TRUE, named distribution representation tracks the
joint distribution of y_t over all horizons.  if FALSE, we track only the
marginal distributions at each horizon.}

\item{quantile_levels}{numeric vector of quantile levels to use for
forecast_representation = "quantile"}

\item{nsim}{integer number of samples to use for
forecast_representation = "sample"}

\item{seed}{either \code{NULL} or an integer that will be used in a call to
\code{set.seed} before simulating the response vectors.  If set, the value is
saved as the "seed" attribute of the returned value.  The default, \code{NULL},
will not change the random generator state, and return \code{.Random.seed}
as the "seed" attribute}

\item{horizon}{number of time steps forwards to simulate}

\item{...}{other arguments passed on to model-specific predict methods}

\item{object}{a model fit of class "lssm", as returned by fit_lssm}
}
\value{
an nsim by horizon matrix with simulated values
}
\description{
This function handles any transformations and differencing that were done
before calling the model-specific predict method.
}
