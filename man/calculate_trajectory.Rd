% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_perturbation.R
\name{calculate_trajectory}
\alias{calculate_trajectory}
\title{Calculate cell trajectory in embedding space}
\usage{
calculate_trajectory(embedding, delta_expr, n_neighbors = 15)
}
\arguments{
\item{embedding}{Original embedding matrix}

\item{delta_expr}{Expression changes}

\item{n_neighbors}{Number of neighbors for trajectory}
}
\value{
Updated embedding with trajectory information
}
\description{
Calculate cell trajectory in embedding space
}
