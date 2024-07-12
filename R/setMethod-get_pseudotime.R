#' @param meta_data meta data
#' @param embeddings embeddings information
#' @param cluster_by choose column used for `slingshot`
#' @param cores CPU cores used for umap
#' @param seed random seed for umap
#' @param start_cluster start_cluster
#' @param end_cluster end_cluster
#' @param verbose verbose
#'
#' @return A list with a matrix and new meta data
#' @export
#'
#' @method get_pseudotime default
#'
#' @rdname get_pseudotime
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' result <- get_pseudotime(example_matrix)
#' }
get_pseudotime.default <- function(
    object,
    meta_data = NULL,
    embeddings = NULL,
    cluster_by = "cluster",
    cores = 1,
    seed = 1,
    start_cluster = NULL,
    end_cluster = NULL,
    verbose = TRUE,
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

  log_message("Running `slingshot`", verbose = verbose)
  slingshot_res <- slingshot::slingshot(
    embeddings,
    clusterLabels = meta_data$cluster,
    start.clus = start_cluster,
    end.clus = end_cluster
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
#' @method get_pseudotime Seurat
#'
#' @rdname get_pseudotime
get_pseudotime.Seurat <- function(
    object,
    assay = "RNA",
    cluster_by = "cluster",
    slot = "data",
    start_cluster = NULL,
    end_cluster = NULL,
    ...) {
  data <- Seurat::GetAssay(object, assay = assay)[slot]
  embeddings <- Seurat::Embeddings(object, reduction = "umap")
  meta_data <- object@meta.data

  result <- get_pseudotime(
    data,
    meta_data = meta_data,
    embeddings = embeddings,
    cluster_by = cluster_by,
    start_cluster = start_cluster,
    end_cluster = end_cluster,
    ...
  )
  meta_data <- result$meta_data[colnames(object), ]
  object@meta.data <- meta_data

  return(object)
}
