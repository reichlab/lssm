% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caret.R
\name{predict_lssm_caret_wrapper}
\alias{predict_lssm_caret_wrapper}
\title{Wrapper around predict.lssm for use with caret::train}
\usage{
predict_lssm_caret_wrapper(modelFit, newdata, ...)
}
\arguments{
\item{newdata}{data from which to generate predictions}

\item{...}{other arguments are ignored}

\item{model_fit}{model fit object of class lssm}
}
\value{
numeric vector of predictive medians with attributes:
\itemize{
\item family is a string with the parametric family, e.g. "norm"
\item other attributes are names of parameters for the parametric family
}
}
\description{
Wrapper around predict.lssm for use with caret::train
}
