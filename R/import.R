#' @import Matrix ggplot2 ggraph tidygraph
#' @importFrom stats family gaussian na.pass
#' @importFrom methods new
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
    "i",
    "x",
    "y",
    "xend",
    "yend",
    "regulator",
    "target",
    "weight",
    "Interaction",
    "name",
    "g",
    "gene",
    "degree",
    "edges",
    "curvetype",
    "betweenness",
    "centrality",
    "communities",
    "weight_new",
    "UMAP_1",
    "UMAP_2",
    "type",
    "n_regions",
    "n_tfs",
    "tf",
    "n_genes",
    "rsq",
    "nvariables",
    "nice",
    "window",
    "EnsDb.Hsapiens.v93.annot.UCSC.hg38", "GRanges", "algorithm",
    "cluster", "comb_dir", "corr", "cutreeDynamicTree", "end_node", "estimate", "findTop",
    "from", "from_node", "gain", "gene_per_tf", "mean_estimate", "mean_padj", "nodes", "normal",
    "padj", "pam", "path", "path_regions", "peak_per_gene", "penalty", "pseudotime", "pval",
    "region", "region_", "seqnames", "start_node", "strand tf_", "tf_per_gene", "to",
    "to_node", "utils_myDist", "val", "verbose", "strand", "tf_"
  )
)
