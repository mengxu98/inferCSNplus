% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{fit_bagging_ridge}
\alias{fit_bagging_ridge}
\title{Fit a bagging ridge regression model as implemented in scikit-learn (python)}
\usage{
fit_bagging_ridge(
  formula,
  data,
  alpha = 1,
  solver = "auto",
  bagging_number = 200L,
  n_jobs = 1,
  p_method = c("wilcox", "t"),
  ...
)
}
\arguments{
\item{formula}{An object of class \code{formula} with a symbolic description
of the model to be fitted.}

\item{data}{A \code{data.frame} containing the variables in the model.}

\item{alpha}{Positive float indicating the regularization strength.}

\item{solver}{Solver to use in the computational routines.
Options include ‘auto’, ‘svd’, ‘cholesky’, ‘lsqr’, ‘sparse_cg’, ‘sag’, ‘saga’.}

\item{bagging_number}{The number of ridge regression model in the bagging.}

\item{n_jobs}{The number of cores used to fit the model.}

\item{p_method}{The test used to calculate p-values. Options are 't' for \code{t.test}, and 'wilcox' for \code{wilcox.test}}

\item{...}{Other parameters for the model fitting function.}
}
\value{
A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
\code{coefs} contains the fitted coefficients.
}
\description{
Fit a bagging ridge regression model as implemented in scikit-learn (python)
}
