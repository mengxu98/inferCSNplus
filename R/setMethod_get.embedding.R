#' @param cores Using n threads for \code{\link[uwot]{umap}}.
#'
#' @export
#'
#' @method get.embedding default
#'
#' @rdname get.embedding
get.embedding.default <- function(
    object,
    cores = 1,
    ...) {
  pca_res <- gmodels::fast.prcomp(object)
  pca_res <- pca_res$x
  rownames(pca_res) <- rownames(object)
  umap_res <- suppressMessages(uwot::umap(object, n_threads = cores))
  rownames(umap_res) <- rownames(object)
  embedding <- list(
    "PCA" = pca_res,
    "UMAP" = umap_res
  )

  return(embedding)
}

#' @export
#'
#' @method get.embedding Seurat
#'
#' @rdname get.embedding
get.embedding.Seurat <- function(object, ...) {
  pca_res <- object@reductions$pca@cell.embeddings
  rownames(pca_res) <- colnames(object)
  umap_res <- object@reductions$umap@cell.embeddings
  rownames(umap_res) <- colnames(object)
  embedding <- list(
    "PCA" = pca_res,
    "UMAP" = umap_res
  )

  return(embedding)
}
