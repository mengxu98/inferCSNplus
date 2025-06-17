#' @include setClass.R
#' @include setGenerics.R

#' @param pfm A \code{PFMatrixList} object with position weight matrices.
#' @param genome A \code{BSgenome} object with the genome of interest.
#' @param motif_tfs A data frame matching motifs with TFs.
#' The first column is assumed to be the name of the motif,
#' the second the name of the TF.
#' @param verbose Display messages.
#'
#' @return A CSNObject object with updated motif info.
#'
#' @rdname find_motifs
#' @export
#' @method find_motifs CSNObject
setMethod(
  f = "find_motifs",
  signature = "CSNObject",
  definition = function(
      object,
      pfm,
      genome,
      motif_tfs = NULL,
      verbose = TRUE,
      ...) {
    params <- Params(object)
    if (is.null(params$peak_assay)) {
      thisutils::log_message(
        "Skipping motif finding as no peak assay was found.",
        verbose = verbose,
        message_type = "warning"
      )
      return(object)
    }

    thisutils::log_message(
      "Adding TF information",
      verbose = verbose
    )
    if (!is.null(motif_tfs)) {
      motif2tf <- motif_tfs
    } else {
      utils::data(motif2tf, envir = environment())
    }

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

    if (length(tfs_use) == 0) {
      stop(
        "None of the provided TFs were found in the dataset.
         Consider providing a custom motif-to-TF map as 'motif_tfs'"
      )
    }
    object@metadata$tfs <- tfs_use
    object@regions@motifs2tfs <- motif2tf[, tfs_use]

    celltypes <- get_attribute(
      object,
      attribute = "celltypes"
    )

    if (params$filter_mode == "variable" && params$filter_by == "aggregate") {
      thisutils::log_message(
        "Processing motifs for all celltypes",
        verbose = verbose
      )
      result <- .process_celltype_motifs(
        object,
        celltypes[1],
        genome,
        pfm,
        verbose = FALSE
      )
      if (!is.null(result)) {
        object@regions@motifs <- purrr::map(
          celltypes,
          ~result
        ) |>
          purrr::set_names(celltypes)
      }
    } else {
      object@regions@motifs <- purrr::map(
        celltypes,
        function(x) {
          .process_celltype_motifs(
            object,
            x,
            genome,
            pfm,
            verbose
          )
        }
      ) |>
        purrr::set_names(celltypes)
    }

    return(object)
  }
)
.process_celltype_motifs <- function(
    object,
    celltype,
    genome,
    pfm,
    verbose) {
  thisutils::log_message(
    "Processing motifs for ", celltype,
    verbose = verbose
  )

  celltype_peaks <- get_attribute(
    object,
    celltypes = celltype,
    attribute = "peaks"
  )

  if (length(celltype_peaks) == 0) {
    thisutils::log_message(
      "no significant peaks found for ", celltype,
      verbose = verbose,
      message_type = "warning"
    )
    return(NULL)
  }

  peak_ranges <- Signac::StringToGRanges(celltype_peaks)
  motif_pos <- suppressWarnings(
    Signac::AddMotifs(
      object = peak_ranges,
      genome = genome,
      pfm = pfm,
      verbose = FALSE
    )
  )

  list(
    motifs = motif_pos,
    n_motifs = ncol(motif_pos),
    peaks = rownames(motif_pos),
    n_peaks = nrow(motif_pos),
    genome = class(genome)[1]
  )
}
