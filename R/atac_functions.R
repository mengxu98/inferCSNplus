#' @param filter_value The value to filter peaks
#'
#' @return
#' @export
#'
#' @rdname peaks.filter
#'
peaks.filter.default <- function(object,
                                 filter_value = 0.9) {
  peaks_table <- Signac::FindTopFeatures(object)
  peaks_table <- peaks_table[peaks_table$percentile > filter_value, ]
  object <- object[rownames(peaks_table), ]
  return(object)
}

#' @param filter_value The value to filter peaks
#'
#' @return
#' @export
#'
#' @rdname peaks.filter matrix
#'
peaks.filter.matrix <- function(object,
                                 filter_value = 0.9) {
  peaks_table <- Signac::FindTopFeatures(object)
  peaks_table <- peaks_table[peaks_table$percentile > filter_value, ]
  object <- object[rownames(peaks_table), ]
  return(object)
}

#' @param filter_value The value to filter peaks
#'
#' @return
#' @export
#'
#' @rdname peaks.filter
#' @method peaks.filter Seurat
#'
peaks.filter.Seurat <- function(object,
                                filter_value = 0.9) {
  DefaultAssay(object) <- "ATAC"
  peaks_table <- Signac::FindTopFeatures(object@assays$ATAC$counts)
  peaks_table <- peaks_table[peaks_table$percentile > filter_value, ]
  object@assays$ATAC@counts <- object@assays$ATAC@counts[rownames(peaks_table), ]
  return(object)
}
