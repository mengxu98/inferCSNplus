#' @export
#'
#' @method get.embedding default
#'
#' @rdname get.embedding
get.embedding.default <- function(object, ...) {
  pca_res <- stats::prcomp(object)
  pca_res <- pca_res$x
  rownames(pca_res) <- rownames(object)
  umap_res <- umap::umap(object)
  rownames(umap_res) <- rownames(object)
  umap_res <- cbind(UMAP1 = umap_res$layout[, 1], UMAP2 = umap_res$layout[, 2])
  embedding <- list(pca_res, umap_res)
  names(embedding) <- c("PCA", "UMAP")

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
  embedding <- list(pca_res, umap_res)
  names(embedding) <- c("PCA", "UMAP")

  return(embedding)
}
