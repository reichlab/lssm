% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossvalidation.R
\name{sarima_param_grid}
\alias{sarima_param_grid}
\title{Define grid of tuning parameters for the arma model}
\usage{
sarima_param_grid(
  y,
  len = NULL,
  search = "grid",
  transformation = "box-cox",
  transform_offset = if (any(x <= 0)) {     -1 * min(x) + 0.49 } else {     0 },
  max_d = 2,
  max_D = 1,
  include_intercept = c(FALSE, TRUE),
  max_p_ar = 5,
  max_q_ma = 5,
  max_P_ar = 2,
  max_Q_ma = 2,
  min_order = 1,
  max_order = 5,
  stationary = 1L
)
}
\arguments{
\item{y}{}

\item{len}{}

\item{search}{}

\item{transform_offset}{}

\item{max_d}{}

\item{max_D}{}

\item{max_p_ar}{}

\item{max_q_ma}{}

\item{max_order}{}

\item{stationary}{}

\item{transform}{}
}
\value{
data frame with columns model, transform, transform_offset, d, D,
p_ar, q_ma, P_ar, Q_ma
}
\description{
Define grid of tuning parameters for the arma model
}
