% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{fit_bayesian_ridge}
\alias{fit_bayesian_ridge}
\title{Fit a Bayesian ridge regression model via scikit-learn}
\usage{
fit_bayesian_ridge(formula, data, ...)
}
\arguments{
\item{formula}{An object of class \emph{\code{formula}} with a symbolic description
of the model to be fitted.}

\item{data}{A \emph{\code{data.frame}} containing the variables in the model.}

\item{...}{Additional parameters passed to scikit-learn's BayesianRidge.}
}
\value{
A list containing two data frames:
\itemize{
\item \emph{\code{gof}} - Model performance metrics
\item \emph{\code{coefs}} - Coefficient estimates with uncertainty measures
}
}
\description{
Fits a Bayesian ridge regression model using scikit-learn's implementation.
}
