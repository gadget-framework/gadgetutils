% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_g3_attributes.R
\name{add_g3_attributes}
\alias{add_g3_attributes}
\title{Add gadget3 attributes to an existing data frame}
\usage{
add_g3_attributes(x, params, warn = 1)
}
\arguments{
\item{x}{A data frame, tibble or data table.}

\item{params}{A list of parameters defining the grouping. See details in \code{mfdb_sample_count}.}

\item{warn}{Integer value to set the handling of warning messages. See the \code{warn} argument in \link[base]{options} for details. In brief, set to negative to ignore warnings, to 1 to print warnings as they occur and to 2 to stop executing and return an error instead of a warning. The last is meant to prevent passing information that may cause unintended behavior in gadget3 models and is recommended.}
}
\value{
Returns \code{x} with g3 attributes and missing columns as instructed by the \code{params} argument.
}
\description{
Adds gadget3 attributes and missing columns to a a data frame ready to be passed to gadget3. Useful for importing data manually to gadget without the mfdb step.
}
