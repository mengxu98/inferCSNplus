% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setGenerics.R, R/find_motifs.R
\name{find_motifs}
\alias{find_motifs}
\alias{find_motifs,CSNObject-method}
\title{Scan for motifs in candidate regions}
\usage{
find_motifs(object, ...)

\S4method{find_motifs}{CSNObject}(object, pfm, genome, motif_tfs = NULL, verbose = TRUE, ...)
}
\arguments{
\item{object}{The input data, a csn object.}

\item{...}{Arguments for other methods}

\item{pfm}{A \code{PFMatrixList} object with position weight matrices.}

\item{genome}{A \code{BSgenome} object with the genome of interest.}

\item{motif_tfs}{A data frame matching motifs with TFs.
The first column is assumed to be the name of the motif,
the second the name of the TF.}

\item{verbose}{Display messages.}
}
\value{
A CSNObject object with updated motif info.
}
\description{
Scan for motifs in candidate regions
}
