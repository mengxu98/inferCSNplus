#' @param pfm A \code{PFMatrixList} object with position weight matrices.
#' @param genome A \code{BSgenome} object with the genome of interest.
#' @param motif_tfs A data frame matching motifs with TFs. The first column is assumed
#' to be the name of the motif, the second the name of the TF.
#' @param regulators regulators.
#' @param verbose Display messages.
#'
#' @return A CSNObject object with updated motif info.
#'
#' @rdname find_motifs
#' @export
#' @method find_motifs CSNObject
find_motifs.CSNObject <- function(
    object,
    pfm,
    genome,
    motif_tfs = NULL,
    regulators = NULL,
    verbose = TRUE,
    ...) {
  params <- Params(object)

  # Add TF info for motifs
  log_message("Adding TF information", verbose = verbose)
  if (!is.null(motif_tfs)) {
    motif2tf <- motif_tfs
  } else {
    utils::data(motif2tf, envir = environment())
  }

  # Spread data frame to sparse matrix
  motif2tf <- motif2tf |>
    dplyr::select("motif" = 1, "tf" = 2) |>
    dplyr::distinct() |>
    dplyr::mutate(val = 1) |>
    tidyr::pivot_wider(
      names_from = "tf",
      values_from = val,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("motif") |>
    as.matrix() |>
    Matrix::Matrix(sparse = TRUE)

  tfs_use <- intersect(
    rownames(GetAssay(object, params$rna_assay)),
    colnames(motif2tf)
  )
  if (!is.null(regulators)) {
    tfs_use <- intersect(
      regulators,
      tfs_use
    )
  }

  if (length(tfs_use) == 0) {
    stop("None of the provided TFs were found in the dataset.
         Consider providing a custom motif-to-TF map as `motif_tfs`")
  }
  object@csn@regions@tfs <- motif2tf[, tfs_use]

  # Find motif positions with Signac/motifmatchr
  cand_ranges <- object@csn@regions@ranges
  motif_pos <- Signac::AddMotifs(
    object = cand_ranges,
    genome = genome,
    pfm = pfm,
    verbose = verbose
  )
  object@csn@regions@motifs <- motif_pos

  return(object)
}
