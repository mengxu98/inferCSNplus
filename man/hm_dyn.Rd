% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_heatmap.R
\name{hm_dyn}
\alias{hm_dyn}
\title{plots results of findDynGenes}
\usage{
hm_dyn(
  matrix,
  dynRes,
  cluster = "cluster",
  topX = 25,
  cRow = FALSE,
  cCol = FALSE,
  limits = c(0, 10),
  toScale = FALSE,
  fontsize_row = 4,
  geneAnn = FALSE,
  anno_colors = NULL,
  show_rownames = TRUE,
  filename = NA,
  width = NA,
  height = NA
)
}
\arguments{
\item{matrix}{expression matrix}

\item{dynRes}{result of running findDynGenes}

\item{cluster}{cluster}

\item{topX}{topX}

\item{cRow}{cRow}

\item{cCol}{cCol}

\item{limits}{limits}

\item{toScale}{toScale}

\item{fontsize_row}{fontsize_row}

\item{geneAnn}{geneAnn}

\item{anno_colors}{anno_colors}

\item{show_rownames}{show_rownames}

\item{filename}{filename}

\item{width}{width}

\item{height}{height}
}
\value{
heatmap list
}
\description{
plots results of findDynGenes
}
