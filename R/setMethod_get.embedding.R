#' @param cores Using n threads for \code{\link[uwot]{umap}}.
#' @param seed Set a random seed, default sets to 1.
#'
#' @export
#'
#' @method get.dimensional default
#'
#' @rdname get.dimensional
get.dimensional.default <- function(
    object,
    cores = 1,
    seed = 1,
    ...) {
  set.seed(seed = seed)

  pca_res <- stats::prcomp(object)
  cell_embeddings <- pca_res$x
  rownames(cell_embeddings) <- rownames(object)
  umap_res <- suppressMessages(
    uwot::umap(object, n_threads = cores)
  )
  rownames(umap_res) <- rownames(object)
  dimensional <- list(
    "PCA" = cell_embeddings,
    "UMAP" = umap_res,
    "feature_loadings" = pca_res$rotation,
    "sdev" = pca_res$sdev
  )
  class(dimensional) <- "VERTOR_dimension"

  return(dimensional)
}


#' @export
#'
#' @method get.dimensional Seurat
#'
#' @rdname get.dimensional
get.dimensional.Seurat <- function(object, ...) {
  dimensional <- list(
    "PCA" = object@reductions$pca@cell.embeddings,
    "UMAP" = object@reductions$umap@cell.embeddings,
    "feature_loadings" = object@reductions$pca@feature.loadings,
    "sdev" = object@reductions$pca@stdev
  )
  class(dimensional) <- "VERTOR_dimension"

  return(dimensional)
}
