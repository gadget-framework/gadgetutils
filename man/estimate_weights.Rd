% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g3_iterative_shortcut.R
\name{estimate_weights}
\alias{estimate_weights}
\title{Shortcut weights}
\usage{
estimate_weights(model, params)
}
\arguments{
\item{model}{A gadget3 model of class 'g3_cpp' or 'g3_r'}

\item{params}{A gadget3 parameter dataframe}
}
\description{
Estimates the weights for gadget3 likelihood components using shortcut (ie not iterative re-weighting) methods
}
\details{
This function attempts to estimate the weights of gadget3 likelihood components. The shortcut weights are calculated by taking the inverse of the residual variance. How the residual variance is calculated depends on the type of likelihood distribution, i.e., whether it is a catch or abundance distribution. The residual variance calculations are performed by the helper functions 'estimate_weights_cdist' (catch distributions) and 'estimate_weights_adist' (abundance distributions)
}
