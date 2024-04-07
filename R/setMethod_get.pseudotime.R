#' @param meta_data meta data
#' @param embeddings embeddings information
#' @param cluster_by choose column used for `slingshot`
#' @param cores CPU cores used for umap
#' @param seed random seed for umap
#'
#' @return A list with a matrix and new meta data
#' @export
#'
#' @method get.pseudotime default
#'
#' @rdname get.pseudotime
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' result <- get.pseudotime(example_matrix)
get.pseudotime.default <- function(
    object,
    meta_data = NULL,
    embeddings = NULL,
    cluster_by = "cluster",
    cores = 1,
    seed = 1,
    ...) {
  if (is.null(meta_data)) {
    cells <- rownames(object)
    if (is.null(cells)) {
      cells <- paste0("cell_", rep(1:nrow(object)))
      rownames(object) <- cells
    }
    meta_data <- data.frame(
      cells = cells,
      cluster = "cluster"
    )
  } else {
    if (!cluster_by %in% colnames(meta_data)) {
      meta_data$cluster <- "cluster"
    } else {
      meta_data$cluster <- meta_data[, cluster_by]
    }
  }

  if (is.null(embeddings)) {
    set.seed(seed)
    embeddings <- suppressMessages(
      uwot::umap(object, n_threads = cores)
    )
  }

  slingshot_res <- slingshot::slingshot(
    embeddings,
    clusterLabels = meta_data$cluster
  )
  pseudotime_res <- slingshot::slingPseudotime(slingshot_res)
  pseudotime_res <- apply(
    pseudotime_res, 2, function(x) {
      normalization(
        x,
        method = "max"
      )
    }
  )
  pseudotime_res <- as.data.frame(pseudotime_res)
  meta_data <- cbind(meta_data, pseudotime_res)
  result <- list(
    matrix = object,
    meta_data = meta_data
  )

  return(result)
}

#' @param assay assay
#' @param slot slot
#'
#' @return Seurat object
#' @export
#'
#' @method get.pseudotime Seurat
#'
#' @rdname get.pseudotime
get.pseudotime.Seurat <- function(
    object,
    assay = "RNA",
    cluster_by = "cluster",
    slot = "data",
    ...) {
  matrix <- Seurat::GetAssay(object, assay = assay)[slot]
  embeddings <- Seurat::Embeddings(object, reduction = "umap")
  meta_data <- object@meta.data

  result <- get.pseudotime(
    matrix,
    meta_data = meta_data,
    embeddings = embeddings,
    cluster_by = cluster_by,
    ...
  )
  meta_data <- result$meta_data[colnames(object), ]
  object@meta.data <- meta_data

  return(object)
}
