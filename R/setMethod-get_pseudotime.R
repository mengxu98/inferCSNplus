#' @param meta_data Input meta data.
#' @param embeddings Embeddings information
#' @param cluster_column The column used for \code{\link[slingshot]{slingshot}}.
#' @param start_cluster The start cluster.
#' @param end_cluster The end cluster.
#' @inheritParams inferCSN
#'
#' @export
#'
#' @method get_pseudotime default
#'
#' @rdname get_pseudotime
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' head(get_pseudotime(example_matrix))
#' }
get_pseudotime.default <- function(
    object,
    meta_data = NULL,
    embeddings = NULL,
    cluster_column = "cluster",
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
    rownames(meta_data) <- cells
  } else {
    if (!cluster_column %in% colnames(meta_data)) {
      meta_data$cluster <- "cluster"
    } else {
      meta_data$cluster <- meta_data[, cluster_column]
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
        x
      )
    }
  )
  pseudotime_res <- as.data.frame(pseudotime_res)
  colnames(pseudotime_res) <- paste0("pseudotime_slingshot", 1:ncol(pseudotime_res))
  pseudotime_res <- cbind.data.frame(cluster = meta_data$cluster, pseudotime_res)

  return(pseudotime_res)
}

#' @param assay The assay used for \code{\link[slingshot]{slingshot}}.
#' @param layer The layer used for \code{\link[slingshot]{slingshot}}.
#' @param reduction The reduction used for \code{\link[slingshot]{slingshot}}.
#'
#' @export
#'
#' @method get_pseudotime Seurat
#'
#' @rdname get_pseudotime
get_pseudotime.Seurat <- function(
    object,
    assay = "RNA",
    layer = "data",
    cluster_column = "cluster",
    reduction = "umap",
    start_cluster = NULL,
    end_cluster = NULL,
    ...) {
  embeddings <- Seurat::Embeddings(object, reduction = reduction)

  result <- get_pseudotime(
    Seurat::GetAssay(object, layer = layer),
    meta_data = object@meta.data,
    embeddings = embeddings,
    cluster_column = cluster_column,
    start_cluster = start_cluster,
    end_cluster = end_cluster,
    ...
  )
  result <- result[colnames(object), ]
  object <- Seurat::AddMetaData(object, result)

  return(object)
}
