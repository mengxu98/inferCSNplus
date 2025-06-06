% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{create_seurat_object}
\alias{create_seurat_object}
\title{Create a Seurat object from simulated data}
\usage{
create_seurat_object(
  sim_data,
  min_cells = 3,
  min_features = 0,
  variable_features = 2000
)
}
\arguments{
\item{sim_data}{Output from \code{\link{simulate_csn_data}}, or a list with
\code{expression_matrix}, \code{peak_matrix}, \code{tfs}, \code{targets},
\code{cell_metadata}, and \code{peak_annotations}.}

\item{min_cells}{Include genes/peaks detected in at least this many cells}

\item{min_features}{Include cells where at least this many genes are detected}

\item{variable_features}{Number of variable features to select}
}
\value{
A Seurat object with RNA and ATAC assays
}
\description{
Create a Seurat object from simulated data
}
\examples{
\dontrun{
sim_data <- simulate_csn_data()
seurat_obj <- create_seurat_object(sim_data)
}
}
