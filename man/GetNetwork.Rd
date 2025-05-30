% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setGenerics.R, R/setMethod.R
\name{GetNetwork}
\alias{GetNetwork}
\alias{GetNetwork,CSNObject-method}
\title{Get network}
\usage{
GetNetwork(object, ...)

\S4method{GetNetwork}{CSNObject}(object, network = DefaultNetwork(object), celltypes = NULL, ...)
}
\arguments{
\item{object}{The input data, a csn object.}

\item{...}{Arguments for other methods}

\item{network}{network}

\item{celltypes}{cell types to analyze, NULL for all cell types}
}
\description{
Get network
}
