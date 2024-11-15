% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setGenerics.R, R/setClass.R
\name{NetworkGraph}
\alias{NetworkGraph}
\alias{NetworkGraph.CSNObject}
\alias{NetworkGraph.RegulatoryNetwork}
\alias{NetworkGraph.Network}
\title{Get network parameters}
\usage{
NetworkGraph(object, ...)

\method{NetworkGraph}{CSNObject}(
  object,
  network = DefaultNetwork(object),
  graph = "module_graph",
  ...
)

\method{NetworkGraph}{RegulatoryNetwork}(
  object,
  network = DefaultNetwork(object),
  graph = "module_graph",
  ...
)

\method{NetworkGraph}{Network}(object, graph = "module_graph", ...)
}
\arguments{
\item{object}{The input data, a csn object.}

\item{...}{Arguments for other methods}

\item{network}{network}

\item{graph}{graph}
}
\description{
Get network parameters
}
