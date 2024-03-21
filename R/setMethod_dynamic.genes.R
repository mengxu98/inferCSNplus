#' Select dynamic genes
#'
#' @param object Expression object.
#' @param pseudotime Pseudotime
#' @param fdr_threshold Threshold of fdr.
#'
#' @return Gene list
#' @export
#'
#' @method dynamic.genes default
#'
#' @rdname dynamic.genes

#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' vector_result <- inferVECTOR(example_matrix)
#' dynamic.genes(
#'   object = t(vector_result@matrix),
#'   pseudotime = vector_result@pseudotime[, 2]
#' )
dynamic.genes.default <- function(
    object,
    pseudotime = NULL,
    fdr_threshold = 0.05,
    ...) {
  if (is.null(pseudotime)) {
    return(rownames(object))
  }
  sorted_genes <- names(
    sort(apply(object, 1, var), decreasing = TRUE)
  )
  object <- object[sorted_genes, ]

  # Fit a GAM model with a loess term for pseudotime
  res <- purrr::map_dfr(
    sorted_genes,
    function(x) {
      data <- data.frame(
        exp = object[x, ],
        t = pseudotime
      )

      gam_model <- suppressWarnings(
        gam::gam(exp ~ gam::lo(t), data = data)
      )
      data.frame(
        gene = x,
        P_value = summary(gam_model)[4][[1]][1, 5]
      )
    }
  )
  res$fdr <- p.adjust(res$P_value, method = "BH")
  genes <- res[res$fdr < fdr_threshold, 1]
  if (length(genes) < 3) {
    genes <- sorted_genes
  }
  return(genes)
}

#' @export
#'
#' @method dynamic.genes VECTOR
#'
#' @rdname dynamic.genes
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' vector_result <- inferVECTOR(example_matrix)
#' dynamic.genes(vector_result)
dynamic.genes.VECTOR <- function(
    object,
    fdr_threshold = 0.05,
    ...) {
  genes <- dynamic.genes(
    t(object@matrix),
    pseudotime = object@pseudotime[, 2],
    fdr_threshold = fdr_threshold
  )

  return(genes)
}

#' @param cluster cluster
#' @param group_by idents
#' @param assay assay used
#'
#' @export
#'
#' @method dynamic.genes Seurat
#'
#' @rdname dynamic.genes
dynamic.genes.Seurat <- function(
    object,
    fdr_threshold = 0.05,
    group_by = NULL,
    cluster = NULL,
    assay = "counts",
    ...) {
  if (!is.null(group_by)) {
    Seurat::Idents(object) <- group_by
  }
  if (!is.null(cluster)) {
    object$cluster <- Seurat::Idents(object)
    object <- subset(object, cluster == cluster)
  }

  genes <- dynamic.genes(
    Matrix::as.matrix(
      switch(
        EXPR = assay,
        "counts" = object@assays$RNA$counts,
        "data" = object@assays$RNA$data
      )
    ),
    pseudotime = object@meta.data$pseudotime,
    fdr_threshold = fdr_threshold
  )

  return(genes)
}
