% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{do_difference}
\alias{do_difference}
\title{Do first-order and seasonal differencing (go from original time series
to differenced time series).}
\usage{
do_difference(y, d = 0, D = 0, frequency = 1)
}
\arguments{
\item{y}{a univariate time series or numeric vector.}

\item{d}{order of first differencing}

\item{D}{order of seasonal differencing}

\item{frequency}{frequency of time series.  Must be provided if y is not
of class "ts" and D > 0.  See the help for stats::ts for more.}
}
\value{
a differenced time series object (of class 'ts'),
padded with leading NAs.
}
\description{
Do first-order and seasonal differencing (go from original time series
to differenced time series).
}
