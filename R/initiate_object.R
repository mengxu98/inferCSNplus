#' @param regions Candidate regions to consider for binding site inference.
#' If \code{NULL}, all peaks regions are considered.
#' @param rna_assay A character vector indicating the name of the gene expression
#' assay in the \code{Seurat} object.
#' @param peak_assay A character vector indicating the name of the chromatin
#' accessibility assay in the \code{Seurat} object.
#' @param exclude_exons Logical. Whether to consider exons for binding site inference.
#'
#' @return A CSNObject object containing a RegulatoryNetwork object.
#'
#' @rdname initiate_object
#' @export
#' @method initiate_object Seurat
initiate_object.Seurat <- function(
    object,
    regions = NULL,
    rna_assay = "RNA",
    peak_assay = "peaks",
    exclude_exons = TRUE,
    ...) {
  gene_annot <- Signac::Annotation(object[[peak_assay]])
  # Through error if NULL
  if (is.null(gene_annot)) {
    stop("Please provide a gene annotation for the ChromatinAssay.")
  }
  peak_ranges <- Signac::StringToGRanges(
    rownames(Seurat::GetAssay(object, assay = peak_assay))
  )

  # Find candidate ranges by intersecting the supplied regions with peaks
  # Per default take all peaks as candidate regions
  if (!is.null(regions)) {
    cand_olaps <- IRanges::findOverlaps(regions, peak_ranges)
    cand_ranges <- IRanges::pintersect(
      peak_ranges[S4Vectors::subjectHits(cand_olaps)],
      regions[S4Vectors::queryHits(cand_olaps)]
    )
  } else {
    cand_ranges <- peak_ranges
  }

  # Exclude exons because they are usually conserved
  if (exclude_exons) {
    exon_ranges <- gene_annot[gene_annot$type == "exon", ]
    names(exon_ranges@ranges) <- NULL
    exon_ranges <- IRanges::intersect(exon_ranges, exon_ranges)
    exon_ranges <- GenomicRanges::GRanges(
      seqnames = exon_ranges@seqnames,
      ranges = exon_ranges@ranges
    )
    cand_ranges <- GenomicRanges::subtract(
      cand_ranges, exon_ranges,
      ignore.strand = TRUE
    ) %>% unlist()
  }

  # Match candidate ranges to peaks
  peak_overlaps <- IRanges::findOverlaps(cand_ranges, peak_ranges)
  peak_matches <- S4Vectors::subjectHits(peak_overlaps)

  regions_obj <- methods::new(
    Class = "Regions",
    ranges = cand_ranges,
    peaks = peak_matches,
    motifs = NULL
  )

  params <- list(
    peak_assay = peak_assay,
    rna_assay = rna_assay,
    exclude_exons = exclude_exons
  )

  csn_obj <- methods::new(
    Class = "RegulatoryNetwork",
    regions = regions_obj,
    params = params
  )

  object <- methods::new(
    Class = "CSNObject",
    csn = csn_obj,
    data = object
  )

  return(object)
}

#' @rdname initiate_object
#' @export
#' @method initiate_object CSNObject
initiate_object.CSNObject <- function(
    object,
    regions = NULL,
    peak_assay = "peaks",
    rna_assay = "RNA",
    exclude_exons = TRUE,
    ...) {
  initiate_object(
    object@data,
    regions = regions,
    peak_assay = peak_assay,
    rna_assay = rna_assay,
    exclude_exons = exclude_exons,
    ...
  )
}

initiate_object_list.CSNObject <- function(
    object,
    cluster = NULL,
    cluster_column = NULL,
    verbose = TRUE,
    peak_assay = "peak",
    pfm = NULL,
    genome = NULL) {
  if (is.null(cluster_column)) {
    return(inferCSN(object)) # ?
  } else {
    clusters <- unique(object@data@meta.data[cluster_column][, 1])

    object_list <- purrr::map(
      clusters, .f = function(x) {
        log_message("Running for cluster: ", x, verbose = verbose)

        object_sub <- subset_object(
          object@data,
          cluster_column = cluster_column,
          cluster = x
        )

        object_sub <- methods::new(
          Class = "CSNObject",
          csn = object@csn,
          data = object_sub
        )

        object_sub <- initiate_object(
          object_sub,
          peak_assay = peak_assay
        )
        # Scan candidate regions for TF binding motifs
        object_sub <- find_motifs(
          object_sub,
          pfm = pfm,
          genome = genome
        )

        return(object_sub)
      }
    )
    names(object_list) <- clusters

    object_list <- methods::new(
      Class = "CSNObjectList",
      data = object_list
    )

    return(object_list)
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
