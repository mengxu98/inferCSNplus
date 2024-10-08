% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setGenerics.R, R/initiate_object.R
\name{initiate_object}
\alias{initiate_object}
\alias{initiate_object.Seurat}
\alias{initiate_object.CSNObject}
\title{Initiate the \code{RegulatoryNetwork} object}
\usage{
initiate_object(object, ...)

\method{initiate_object}{Seurat}(
  object,
  regions = NULL,
  rna_assay = "RNA",
  peak_assay = "peaks",
  exclude_exons = TRUE,
  ...
)

\method{initiate_object}{CSNObject}(
  object,
  regions = NULL,
  peak_assay = "peaks",
  rna_assay = "RNA",
  exclude_exons = TRUE,
  ...
)
}
\arguments{
\item{object}{The input data, a seurat object.}

\item{...}{Arguments for other methods}

\item{regions}{Candidate regions to consider for binding site inference.
If \code{NULL}, all peaks regions are considered.}

\item{rna_assay}{A character vector indicating the name of the gene expression
assay in the \code{Seurat} object.}

\item{peak_assay}{A character vector indicating the name of the chromatin
accessibility assay in the \code{Seurat} object.}

\item{exclude_exons}{Logical. Whether to consider exons for binding site inference.}
}
\value{
A CSNObject object containing a RegulatoryNetwork object.
}
\description{
Initiate the \code{RegulatoryNetwork} object
}
