#' @import Matrix ggplot2 ggraph tidygraph
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict cor na.omit sd
#' @importFrom utils methods
#' @importFrom stats family gaussian na.pass
#' @importFrom methods as is new
#' @importClassesFrom Signac Motif
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SeuratObject Seurat
NULL

#' @import sparseMatrixStats
summary_fun <- list(
  "mean" = sparseMatrixStats::colMeans2,
  "median" = sparseMatrixStats::colMedians,
  "max" = sparseMatrixStats::colMaxs,
  "min" = sparseMatrixStats::colMins,
  "count" = sparseMatrixStats::colCounts,
  "any" = sparseMatrixStats::colAnys,
  "all" = sparseMatrixStats::colAlls,
  "sd" = sparseMatrixStats::colSds,
  "mad" = sparseMatrixStats::colMads
)

utils::globalVariables(
  c(
    ".",
    "aggregate",
    "algorithm",
    "betweenness",
    "centrality",
    "celltype",
    "cluster",
    "comb_dir",
    "communities",
    "components",
    "corr",
    "curvetype",
    "cutreeDynamicTree",
    "degree",
    "edges",
    "end_node",
    "EnsDb.Hsapiens.v93.annot.UCSC.hg38",
    "estimate",
    "findTop",
    "from",
    "from_node",
    "g",
    "gain",
    "gene",
    "gene_per_tf",
    "GRanges",
    "i",
    "id",
    "Interaction",
    "label_genes",
    "mean_estimate",
    "mean_padj",
    "n_genes",
    "n_regions",
    "n_tfs",
    "name",
    "nice",
    "nodes",
    "normal",
    "nvariables",
    "padj",
    "pam",
    "path",
    "path_regions",
    "peak_per_gene",
    "penalty",
    "pseudotime",
    "pval",
    "region",
    "region_",
    "regulator",
    "rsq",
    "seqnames",
    "start_node",
    "strand tf_",
    "strand",
    "target",
    "targets_num",
    "tf",
    "tf_",
    "tf_per_gene",
    "to",
    "to_node",
    "type",
    "UMAP_1",
    "UMAP_2",
    "utils_myDist",
    "val",
    "verbose",
    "weight",
    "weight_new",
    "window",
    "x",
    "xend",
    "y",
    "yend"
  )
)
