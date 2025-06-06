% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{fit_xgb}
\alias{fit_xgb}
\title{Fit a gradient boosting regression model with XGBoost}
\usage{
fit_xgb(
  formula,
  data,
  params = list(max_depth = 3, eta = 0.01, objective = "reg:squarederror"),
  nrounds = 1000,
  nthread = -1,
  ...
)
}
\arguments{
\item{formula}{An object of class \emph{\code{formula}} with a symbolic description
of the model to be fitted.}

\item{data}{A \emph{\code{data.frame}} containing the variables in the model.}

\item{params}{A list of XGBoost parameters including:
\itemize{
\item \emph{\code{max_depth}} - Maximum tree depth
\item \emph{\code{eta}} - Learning rate
\item \emph{\code{objective}} - Loss function to optimize
}}

\item{nrounds}{Maximum number of boosting iterations.}

\item{nthread}{Number of parallel threads used (-1 for all cores).}

\item{...}{Additional parameters passed to \emph{\code{xgboost::xgboost}}.}
}
\value{
A list containing two data frames:
\itemize{
\item \emph{\code{metrics}} - Model performance metrics
\item \emph{\code{coefficients}} - Feature importance measures
}
}
\description{
Fits a gradient boosted tree model using the XGBoost implementation.
}
