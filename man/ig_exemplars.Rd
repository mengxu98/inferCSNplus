% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_heatmap.R
\name{ig_exemplars}
\alias{ig_exemplars}
\title{make a graph of the regulators, top targets, selecting only top XX targets each}
\usage{
ig_exemplars(grnTab, geneDF, tfList, topX = 5, posOnly = TRUE)
}
\arguments{
\item{grnTab}{grnTab}

\item{geneDF}{geneDF}

\item{tfList}{tfList}

\item{topX}{topX}

\item{posOnly}{posOnly}
}
\value{
iGraph object
}
\description{
make a graph of the regulators, top targets, selecting only top XX targets each
}
