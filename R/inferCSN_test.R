infer_csn_test.CSNObject <- function(
    object,
    cluster = NULL,
    cluster_column = NULL) {
  if (is.null(cluster_column)) {
    return(inferCSN(object))
  } else {
    clusters <- unique(object@data@meta.data[cluster_column][, 1])
    object_list <- purrr::map(
      clusters, .f = function(x) {
        object_sub <- subset_object(
          object@data,
          cluster_column = cluster_column,
          cluster = x
        )
      }
    )
  }
}

subset_object <- function(
    object,
    cluster = NULL,
    cluster_column = NULL){
  object[, object@meta.data[cluster_column] == cluster]
}



#' Standardize the values by rows
#'
#' This function standardize the values for a matrix by its rows
#'
#' @return a gene by cell (or pseudotime) expression matrix
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param data a matrix.
scale_matrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data, 1, sd)
  (data - cm) / csd
}
