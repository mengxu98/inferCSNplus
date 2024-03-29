% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.visualization.R
\name{network.heatmap}
\alias{network.heatmap}
\title{The heatmap of network}
\usage{
network.heatmap(
  weight_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE,
  show_names = FALSE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  heatmap_color = c("#1966ad", "white", "#bb141a"),
  legend_name = "Weight",
  row_title = "Regulators",
  abs_weight = FALSE
)
}
\arguments{
\item{weight_table}{The weight data table of network.}

\item{regulators}{Regulators list.}

\item{targets}{Targets list.}

\item{switch_matrix}{Logical value, whether to weight data table to matrix.}

\item{show_names}{Logical value, whether to show names of row and column.}

\item{heatmap_size}{The size of heatmap, default set to 5.}

\item{heatmap_height}{The height of heatmap.}

\item{heatmap_width}{The width of heatmap.}

\item{heatmap_title}{The title of heatmap.}

\item{heatmap_color}{Colors of heatmap.}

\item{legend_name}{The name of legend.}

\item{row_title}{The title of row.}

\item{abs_weight}{Logical value, whether to perform absolute value on weights.}
}
\value{
Return a heatmap of ggplot2 object
}
\description{
The heatmap of network
}
\examples{
library(inferCSN)
data("example_matrix")
data("example_ground_truth")
weight_table <- inferCSN(example_matrix)

p1 <- network.heatmap(
  example_ground_truth[, 1:3],
  heatmap_title = "Ground truth",
  legend_name = "Ground truth"
)
p2 <- network.heatmap(
  weight_table,
  heatmap_title = "inferCSN",
  legend_name = "inferCSN"
)
ComplexHeatmap::draw(p1 + p2)

p3 <- network.heatmap(
  weight_table,
  heatmap_title = "inferCSN",
  legend_name = "Weight1",
  heatmap_color = c("#20a485", "#410054", "#fee81f")
)
p4 <- network.heatmap(
  weight_table,
  heatmap_title = "inferCSN",
  legend_name = "Weight2",
  heatmap_color = c("#20a485", "white", "#fee81f")
)
ComplexHeatmap::draw(p3 + p4)

network.heatmap(
  weight_table,
  heatmap_title = "inferCSN",
  show_names = TRUE
)

network.heatmap(
  weight_table,
  regulators = c("g1", "g2"),
  heatmap_title = "inferCSN",
  show_names = TRUE
)

network.heatmap(
  weight_table,
  targets = c("g1", "g2"),
  heatmap_title = "inferCSN",
  show_names = TRUE
)

network.heatmap(
  weight_table,
  regulators = c("g1", "g3", "g5"),
  targets = c("g3", "g6", "g9"),
  heatmap_title = "inferCSN",
  show_names = TRUE
)
}
