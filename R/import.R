#' @import Matrix ggplot2 ggraph patchwork tidygraph ggnetwork
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef predict
#' @importFrom utils methods
#' @importFrom stats family gaussian na.pass
#' @importFrom methods as is new
#' @importClassesFrom Signac Motif
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SeuratObject Seurat
NULL

makeNindexFromArrayViewport <- function(
    viewport,
    expand.RangeNSBS = FALSE) {
  viewport_ranges <- IRanges::ranges(viewport)
  viewport_dim <- dim(viewport)
  viewport_refdim <- S4Arrays::refdim(viewport)
  ndim <- length(viewport_dim)
  Nindex <- vector("list", length = ndim)
  is_not_missing <- viewport_dim < viewport_refdim
  if (expand.RangeNSBS) {
    expand_idx <- which(is_not_missing)
  } else {
    viewport_starts <- stats::start(viewport_ranges)
    viewport_ends <- stats::end(viewport_ranges)
    is_width1 <- viewport_dim == 1L
    expand_idx <- which(is_not_missing & is_width1)
    RangeNSBS_idx <- which(is_not_missing & !is_width1)
    Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx, function(i) {
      range_start <- viewport_starts[[i]]
      range_end <- viewport_ends[[i]]
      upper_bound <- viewport_refdim[[i]]
      new2("RangeNSBS",
        subscript = c(range_start, range_end),
        upper_bound = upper_bound, check = FALSE
      )
    })
  }
  if (length(expand_idx) != 0L) {
    Nindex[expand_idx] <- as.list(as(
      viewport_ranges[expand_idx],
      "CompressedIntegerList"
    ))
  }
  Nindex
}

summary_fun <- list(
  "makeNindexFromArrayViewport" = makeNindexFromArrayViewport
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
    "DR1",
    "DR2",
    "edges",
    "end_node",
    "EnsDb.Hsapiens.v93.annot.UCSC.hg38",
    "Expression",
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
    "Iteration",
    "k",
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
    "new2",
    "nice",
    "nodes",
    "normal",
    "nvariables",
    "P_k",
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
    "r_squared",
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
