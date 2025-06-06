% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamic_network.R
\name{define_epochs}
\alias{define_epochs}
\title{Define epochs}
\usage{
define_epochs(
  dynamic_object,
  matrix,
  method = "pseudotime",
  num_epochs = 2,
  pseudotime_cuts = NULL,
  group_assignments = NULL,
  p_value = 0.05,
  winSize = 2
)
}
\arguments{
\item{dynamic_object}{results of running findDynGenes, or a list of results of running findDynGenes per path. If list, names should match names of matrix.}

\item{matrix}{genes-by-cells expression matrix, or a list of expression matrices per path. If list, names should match names of dynamic_object.}

\item{method}{method to define epochs. Either "pseudotime", "cell_order", "group", "con_similarity", "kmeans", "hierarchical"}

\item{num_epochs}{number of epochs to define. Ignored if epoch_transitions, pseudotime_cuts, or group_assignments are provided.}

\item{pseudotime_cuts}{vector of pseudotime cutoffs. If NULL, cuts are set to max(pseudotime)/num_epochs.}

\item{group_assignments}{a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch}

\item{p_value}{p_value}

\item{winSize}{winSize}
}
\value{
updated list of dynamic_object with epoch column included in dynamic_object$cells
}
\description{
Define epochs
}
