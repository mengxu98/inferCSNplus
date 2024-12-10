#' @import Matrix ggplot2 ggraph patchwork tidygraph ggnetwork
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef predict
#' @importFrom utils methods
#' @importFrom stats family gaussian na.pass
#' @importFrom methods as is new
#' @importFrom DelayedArray makeNindexFromArrayViewport
#' @importClassesFrom Signac Motif
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SeuratObject Seurat
NULL

summary_fun <- list(
  "makeNindexFromArrayViewport" = DelayedArray::makeNindexFromArrayViewport
)

utils::globalVariables(
  c(
    ".",
    "Actual",
    "Category",
    "Metric",
    "P_value",
    "P_value_contrary",
    "Predicted",
    "Type",
    "UMAP_1",
    "UMAP_2",
    "Value",
    "aggregate",
    "algorithm",
    "betweenness",
    "celltype",
    "centrality",
    "cluster",
    "comb_dir",
    "communities",
    "components",
    "corr",
    "count",
    "curvetype",
    "degree",
    "edges",
    "end_node",
    "EnsDb.Hsapiens.v93.annot.UCSC.hg38",
    "coefficient",
    "from",
    "from_node",
    "fold.enrichment",
    "g",
    "gain",
    "gene",
    "gene_per_tf",
    "GRanges",
    "i",
    "id",
    "Interaction",
    "label_genes",
    "lda_data",
    "mean_corr",
    "mean_estimate",
    "mean_padj",
    "mean_weight",
    "motif",
    "motif2tf",
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
    "pvalue",
    "r2",
    "region",
    "region_",
    "regulator",
    "rsq",
    "seqnames",
    "start_node",
    "strand tf_",
    "strand",
    "sum_weight",
    "supercell_GE",
    "target",
    "targets_num",
    "tf",
    "tf_",
    "tf_per_gene",
    "to",
    "to_node",
    "type",
    "v",
    "v_new",
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
