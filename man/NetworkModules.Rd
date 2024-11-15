% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setGenerics.R, R/setClass.R
\name{NetworkModules}
\alias{NetworkModules}
\alias{NetworkModules.CSNObject}
\alias{NetworkModules.RegulatoryNetwork}
\alias{NetworkModules.Network}
\title{Get TF modules}
\usage{
NetworkModules(object, ...)

\method{NetworkModules}{CSNObject}(object, network = DefaultNetwork(object), ...)

\method{NetworkModules}{RegulatoryNetwork}(object, network = DefaultNetwork(object), ...)

\method{NetworkModules}{Network}(object, ...)
}
\arguments{
\item{object}{The input data, a csn object.}

\item{...}{Arguments for other methods}

\item{network}{network}
}
\description{
Get TF modules
}
