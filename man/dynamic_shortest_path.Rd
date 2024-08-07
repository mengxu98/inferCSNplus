% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_analysis.R
\name{dynamic_shortest_path}
\alias{dynamic_shortest_path}
\title{Function to return shortest path from 1 regulator to 1 target in a dynamic network}
\usage{
dynamic_shortest_path(
  network_table,
  regulator,
  target,
  weight_column = "weight",
  compare_to_average = FALSE
)
}
\arguments{
\item{network_table}{a dyanmic network}

\item{regulator}{the starting regulator}

\item{target}{the end regulator}

\item{weight_column}{column name in network_table with edge weights that will be converted to distances}

\item{compare_to_average}{if TRUE will compute normalized against average path length}
}
\value{
shortest path, distance, normalized distance, and action
}
\description{
Function to return shortest path from 1 regulator to 1 target in a dynamic network
}
