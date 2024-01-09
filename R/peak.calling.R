#' peak.calling
#'
#' @param macs2.path path to macs2
#' @param fragments fragments
#' @inheritParams inferCSN
#'
#' @return seurat object
#' @export
peak.calling <- function(
    object,
    fragments = NULL,
    macs2_path = NULL,
    verbose = FALSE) {
  if (is.null(fragments)) stop("Please input fragments!")
  if (is.null(macs2_path)) message("Please give the path to macs2!")
  # https://macs3-project.github.io/MACS/

  object$cluster <- Seurat::Idents(object)
  DefaultAssay(object) <- "ATAC"
  peaks <- Signac::CallPeaks(
    object = object,
    group.by = "cluster",
    macs2.path = macs2_path
  )

  object@assays$ATAC@counts <- Signac::FeatureMatrix(
    fragments = fragments,
    features = peaks
  )

  if (verbose) message("Peak calling finished.")
}
