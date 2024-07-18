#' @param cores Using n threads for \code{\link[uwot]{umap}}.
#' @param seed Set a random seed, default sets to 1.
#' @param dims The number of dimensions to use.
#'
#' @export
#'
#' @method get_embedding default
#'
#' @rdname get_embedding
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' dimensional_information <- get_embedding(example_matrix)
#' }
setMethod(
  f = "get_embedding",
  signature = signature(object = "matrix"),
  definition = function(
      object,
      cores = 1,
      seed = 1,
      dims = 10,
      ...) {
    set.seed(seed = seed)

    pca_res <- stats::prcomp(object)
    cell_embeddings <- pca_res$x
    rownames(cell_embeddings) <- rownames(object)
    umap_res <- suppressMessages(
      uwot::umap(object, n_threads = cores)
    )
    umap_res <- suppressMessages(
      uwot::umap(cell_embeddings[, 1:dims], n_threads = cores)
    )
    rownames(umap_res) <- rownames(object)

    return(
      list(
        "pca" = cell_embeddings,
        "umap" = umap_res,
        "feature_loadings" = pca_res$rotation,
        "sdev" = pca_res$sdev
      )
    )
  }
)

#' @export
#'
#' @method get_embedding Seurat
#' @param reduction The reduction to use.
#' @param dims The number of dimensions to use.
#'
#' @rdname get_embedding
setMethod(
  f = "get_embedding",
  signature = signature(object = "Seurat"),
  definition = function(
      object,
      reduction = "umap",
      dims = 2,
      ...) {
    return(
      list(
        "pca" = Seurat::Embeddings(object, reduction = "pca"),
        "umap" = Seurat::Embeddings(object, reduction = reduction)[, 1:dims],
        "feature_loadings" = object@reductions$pca@feature.loadings,
        "sdev" = object@reductions$pca@stdev
      )
    )
  }
)
