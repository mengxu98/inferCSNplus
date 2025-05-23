% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_regression_model.R
\name{fit_srm}
\alias{fit_srm}
\title{Sparse regression model}
\usage{
fit_srm(
  x,
  y,
  cross_validation = FALSE,
  seed = 1,
  penalty = "L0",
  regulators_num = ncol(x),
  n_folds = 5,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{The matrix of regulators.}

\item{y}{The vector of target.}

\item{cross_validation}{Logical value, default is \emph{\code{FALSE}}, whether to use cross-validation.}

\item{seed}{The random seed for cross-validation, default is \emph{\code{1}}.}

\item{penalty}{The type of regularization, default is \emph{\code{L0}}.
This can take either one of the following choices: \emph{\code{L0}}, \emph{\code{L0L1}}, and \emph{\code{L0L2}}.
For high-dimensional and sparse data, \emph{\code{L0L2}} is more effective.}

\item{regulators_num}{The number of non-zore coefficients, this value will affect the final performance.
The maximum support size at which to terminate the regularization path.}

\item{n_folds}{The number of folds for cross-validation, default is \emph{\code{5}}.}

\item{verbose}{Logical value, default is \emph{\code{TRUE}}, whether to print progress messages.}

\item{...}{Parameters for other methods.}
}
\value{
A list of the sparse regression model.
The list has the following components:
\item{model}{The sparse regression model.}
\item{metrics}{A list of metrics.}
\item{coefficients}{A list of coefficients.}
}
\description{
Sparse regression model
}
\examples{
data("example_matrix")
fit_srm(
  x = example_matrix[, -1],
  y = example_matrix[, 1]
)
}
