% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_perturbation.R
\name{predictPerturbation}
\alias{predictPerturbation}
\alias{predictPerturbation,Network-method}
\alias{predictPerturbation,CSNObject-method}
\title{Predict gene expression under TF perturbation}
\usage{
predictPerturbation(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 5,
  scale_factor = 0.1,
  verbose = TRUE
)

\S4method{predictPerturbation}{Network}(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 3,
  scale_factor = 0.1,
  verbose = TRUE
)

\S4method{predictPerturbation}{CSNObject}(object, perturb_tfs, use_weight = FALSE)
}
\arguments{
\item{object}{A Network or CSNObject object}

\item{perturb_tfs}{Named vector of TF perturbation states}

\item{use_weight}{Logical, whether to use network weights instead of coefficients}

\item{n_iter}{Number of iterations for perturbation simulation}

\item{scale_factor}{Scale factor for perturbation effect (0-1)}

\item{verbose}{Whether to show progress messages}
}
\value{
Updated object with perturbation predictions
}
\description{
Predict gene expression under TF perturbation
}
