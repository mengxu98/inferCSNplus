% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamic_genes.R
\name{assign_network}
\alias{assign_network}
\title{Assigns genes to epochs}
\usage{
assign_network(
  matrix,
  dynamic_object,
  method = "active_expression",
  p_value = 0.05,
  pThresh_DE = 0.05,
  active_thresh = 0.33,
  toScale = FALSE,
  forceGenes = TRUE
)
}
\arguments{
\item{matrix}{genes-by-cells expression matrix}

\item{dynamic_object}{individual path result of running define_epochs}

\item{method}{method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)}

\item{p_value}{pval threshold if gene is dynamically expressed}

\item{pThresh_DE}{pval if gene is differentially expressed. Ignored if method is active_expression.}

\item{active_thresh}{value between 0 and 1. Percent threshold to define activity}

\item{toScale}{whether or not to scale the data}

\item{forceGenes}{whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.}
}
\value{
epochs a list detailing genes active in each epoch
}
\description{
Assigns genes to epochs
}
