#' peaks.filter
#' @param object matrix or seurat object
#' @param filter_value filter_value
#'
#' @export
peaks.filter <- function(
    object,
    filter_value = 0.9) {
  if (is(object, "Seurat")) {
    DefaultAssay(object) <- "ATAC"
    peaks_table <- Signac::FindTopFeatures(object@assays$ATAC$counts)
    peaks_table <- peaks_table[peaks_table$percentile > filter_value, ]
    object@assays$ATAC@counts <- object@assays$ATAC@counts[rownames(peaks_table), ]
  } else {
    peaks_table <- Signac::FindTopFeatures(object)
    peaks_table <- peaks_table[peaks_table$percentile > filter_value, ]
    object <- object[rownames(peaks_table), ]
  }

  return(object)
}
