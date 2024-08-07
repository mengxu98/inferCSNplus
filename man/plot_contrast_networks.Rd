% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_visualization.R
\name{plot_contrast_networks}
\alias{plot_contrast_networks}
\title{plot_contrast_networks}
\usage{
plot_contrast_networks(
  network_table,
  degree_value = 0,
  weight_value = 0,
  legend_position = "bottom"
)
}
\arguments{
\item{network_table}{The weight data table of network.}

\item{degree_value}{degree_value}

\item{weight_value}{weight_value}

\item{legend_position}{The position of legend.}
}
\value{
Return a ggplot2 object
}
\description{
plot_contrast_networks
}
\examples{
data("example_matrix")
network_table <- inferCSN(example_matrix)
plot_contrast_networks(network_table[1:50, ])
}
