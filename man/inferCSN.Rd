% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inferCSN.R, R/inferCSN-methods.R
\docType{methods}
\name{inferCSN}
\alias{inferCSN}
\alias{inferCSN,matrix-method}
\alias{inferCSN,sparseMatrix-method}
\alias{inferCSN,data.frame-method}
\alias{inferCSN,Seurat-method}
\alias{inferCSN,CSNObject-method}
\alias{inferCSN,CSNObjectList-method}
\title{\strong{infer}ring \strong{C}ell-\strong{S}pecific gene regulatory \strong{N}etwork}
\usage{
inferCSN(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

\S4method{inferCSN}{matrix}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

\S4method{inferCSN}{sparseMatrix}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

\S4method{inferCSN}{data.frame}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

\S4method{inferCSN}{Seurat}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = FALSE,
  aggregate = TRUE,
  k_neigh = 50,
  atacbinary = TRUE,
  max_overlap = 0.8,
  reduction_name = NULL,
  size_factor_normalize = FALSE,
  genome_info = NULL,
  high_corr_cutoff = NULL,
  low_corr_cutoff = NULL,
  rescued = TRUE,
  ...
)

\S4method{inferCSN}{CSNObject}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  genes = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  tf_cor = 0.1,
  peak_cor = 0,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("srm", "glm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge",
    "bayesian_ridge"),
  alpha = 0.5,
  family = "gaussian",
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  ...
)

\S4method{inferCSN}{CSNObjectList}(
  object,
  penalty = "L0",
  algorithm = "CD",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  cores = 1,
  verbose = TRUE,
  genes = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  tf_cor = 0.1,
  peak_cor = 0,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("srm", "glm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge",
    "bayesian_ridge"),
  alpha = 0.5,
  family = "gaussian",
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  clusters = NULL,
  ...
)
}
\arguments{
\item{object}{The input data for \emph{\code{inferCSN}}.}

\item{penalty}{The type of regularization, default is \emph{\code{L0}}.
This can take either one of the following choices: \emph{\code{L0}}, \emph{\code{L0L1}}, and \emph{\code{L0L2}}.
For high-dimensional and sparse data, \emph{\code{L0L2}} is more effective.}

\item{algorithm}{The type of algorithm used to minimize the objective function, default is \emph{\code{CD}}.
Currently \emph{\code{CD}} and \emph{\code{CDPSI}} are supported.
The \emph{\code{CDPSI}} algorithm may yield better results, but it also increases running time.}

\item{cross_validation}{Logical value, default is \emph{\code{FALSE}}, whether to use cross-validation.}

\item{seed}{The random seed for cross-validation, default is \emph{\code{1}}.}

\item{n_folds}{The number of folds for cross-validation, default is \emph{\code{10}}.}

\item{subsampling}{The percent of all samples used for \code{\link{sparse_regression}}, default is \emph{\code{1}}.}

\item{r_threshold}{Threshold of \eqn{R^2} or correlation coefficient, default is \emph{\code{0}}.}

\item{regulators}{The regulator genes for which to infer the regulatory network.}

\item{targets}{The target genes for which to infer the regulatory network.}

\item{regulators_num}{The number of non-zore coefficients, this value will affect the final performance.
The maximum support size at which to terminate the regularization path.
Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.}

\item{cores}{The number of cores to use for parallelization with \code{\link[foreach]{foreach}}, default is \emph{\code{1}}.}

\item{verbose}{Logical value, default is \emph{\code{TRUE}}, whether to print progress messages.}

\item{...}{Parameters for other methods.}

\item{aggregate}{Logical, whether to aggregate the data.}

\item{k_neigh}{Number of cells to be aggregated per cluster.}

\item{atacbinary}{Logical, whether the aggregated scATAC-seq data need binary}

\item{max_overlap}{The maximum overlapping ratio of two clusters.}

\item{reduction_name}{The reduction name of extracting the cell coordinates used for aggregating.}

\item{size_factor_normalize}{Logical, should accessibility values be normalized by size factor.}

\item{genome_info}{Genome information.}

\item{high_corr_cutoff}{atacbinary}

\item{low_corr_cutoff}{atacbinary}

\item{rescued}{Logical}

\item{genes}{A character vector with the target genes to consider for CSN inference.
Takes all \code{VariableFeatures} in the object per default.}

\item{network_name}{network_name.}

\item{peak_to_gene_method}{Character specifying the method to
link peak overlapping motif regions to nearby genes. One of \code{Signac} or \code{GREAT}.}

\item{upstream}{Integer defining the distance upstream of the gene to consider as potential regulatory region.}

\item{downstream}{Integer defining the distance downstream of the gene to consider as potential regulatory region.}

\item{extend}{Integer defining the distance from the upstream and downstream of the basal regulatory region.
Only used of `peak_to_gene_method = 'GREAT'`.}

\item{only_tss}{Logical. Measure distance from the TSS (\code{TRUE}) or from the entire gene body (\code{FALSE}).}

\item{peak_to_gene_domains}{\code{GenomicRanges} object with regulatory regions for each gene.}

\item{tf_cor}{Threshold for TF - target gene correlation.}

\item{peak_cor}{Threshold for binding peak - target gene correlation.}

\item{aggregate_rna_col}{aggregate_rna_col}

\item{aggregate_peaks_col}{aggregate_peaks_col}

\item{method}{A character string indicating the method to fit the model.
* \code{'srm'} - Sparse Regression Model.
* \code{'glm'} - Generalized Liner Model with \code{\link[stats]{glm}}.
* \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
* \code{'brms'} - Bayesian Regression Models using \code{\link[brms]{brms-package}}.
* \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
* \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.
* \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.}

\item{alpha}{The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.}

\item{family}{A description of the error distribution and link function to be used in the model.
See \code{\link[stats]{family}} for mode details.}

\item{interaction_term}{The interaction term to use in the model between TF and binding site.
* \code{'+'} for additive interaction.
* \code{':'} for 'multiplicative' interaction.
* \code{'*'} for crossing interaction, i.e. additive AND 'multiplicative'.
For more info, see \code{\link[stats]{formula}}}

\item{adjust_method}{Method for adjusting p-values.}

\item{scale}{Logical. Whether to z-transform the expression and accessibility matrices.}

\item{clusters}{clusters.}
}
\value{
A data table of regulator-target regulatory relationships

A CSNObject.

A CSNObject_list.
}
\description{
\strong{infer}ring \strong{C}ell-\strong{S}pecific gene regulatory \strong{N}etwork
}
\examples{
data("example_matrix")
network_table_1 <- inferCSN(
  example_matrix
)
head(network_table_1)

network_table_2 <- inferCSN(
  example_matrix,
  cores = 2
)

identical(
  network_table_1,
  network_table_2
)

inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g0")
)
\dontrun{
data("example_matrix")
network_table <- inferCSN(example_matrix)
head(network_table)

network_table_sparse_1 <- inferCSN(
  as(example_matrix, "sparseMatrix")
)
head(network_table_sparse_1)

network_table_sparse_2 <- inferCSN(
  as(example_matrix, "sparseMatrix"),
  cores = 2
)
identical(
  network_table,
  network_table_sparse_1
)

identical(
  network_table_sparse_1,
  network_table_sparse_2
)

plot_scatter(
  data.frame(
    network_table$weight,
    network_table_sparse_1$weight
  ),
  legend_position = "none"
)

plot_weight_distribution(
  network_table
) + plot_weight_distribution(
  network_table_sparse_1
)
}
\dontrun{
data("promoter_regions_hg38")
seurat_object <- inferCSN(
  seurat_object,
  genome_info = promoter_regions_hg38
)
}
}
