% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g3_jitter.R
\name{g3_jitter}
\alias{g3_jitter}
\title{Jitter analysis}
\usage{
g3_jitter(
  gd,
  outdir = "JITTER",
  model,
  params,
  njits = 10,
  jitter_fraction = 0.1,
  pattern_to_ignore = NULL,
  within_bounds = TRUE,
  use_parscale = TRUE,
  method = "BFGS",
  control = list(),
  serial_compile = FALSE,
  mc.cores = parallel::detectCores()
)
}
\arguments{
\item{gd}{Directory to store output}

\item{outdir}{Directory name within gd to store run outputs}

\item{model}{A G3 model, produced by g3_to_tmb() or g3_to_r()}

\item{params}{Initial parameters to use with the model, this should be a TMB parameter template i.e. attr(tmb_model, 'parameter_template')}

\item{njits}{Number of jitters to run, defaults to 10}

\item{jitter_fraction}{The fraction of jittering for a value}

\item{pattern_to_ignore}{Regular expression of parameters to avoid jittering}

\item{within_bounds}{Logical, if TRUE, jittered values that fall outside parameter bounds (params$lower, params$upper) will be adjusted to fall within the bounds}

\item{use_parscale}{Logical indicating whether optim(control$parscale) should be used}

\item{method}{The optimisation method, see \code{\link[stats]{optim}}}

\item{control}{List of control options for optim, see \code{\link[stats]{optim}}}

\item{serial_compile}{g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues}

\item{mc.cores}{The number of cores to use, defaults to the number available}
}
\value{
A list of optimised parameter data frames (one for each jitter)
}
\description{
\code{g3_jitter} performs a jitter analysis.
}
