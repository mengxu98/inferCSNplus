% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_comparison.R
\name{JI_across_topregs}
\alias{JI_across_topregs}
\title{Computes Jaccard similarity between top regulators in two sets of networks across a range of top X regulators}
\usage{
JI_across_topregs(
  netlist1,
  netlist2,
  n_regs = 3:15,
  func = "mean",
  method = "pagerank",
  weight_column = "zscore",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
)
}
\arguments{
\item{netlist1}{list of grnDFs}

\item{netlist2}{list of grnDFs}

\item{n_regs}{a vector indicating which values of top regulators to scan across}

\item{func}{func}

\item{method}{method to find top regulators. Currently only supports "pagerank"}

\item{weight_column}{column name in grnDFs containing edge weights}

\item{compare_within_netlist1}{whether or not to do pairwise comparisons between networks in netlist1}

\item{compare_within_netlist2}{whether or not to do pairwise comparisons between networks in netlist2}
}
\value{
dataframe of Jaccard similarities of top regulators
}
\description{
Computes Jaccard similarity between top regulators in two sets of networks across a range of top X regulators
}
