PlotTypeSpecificMotifs <- function(
    object,
    cell_type,
    n_motifs = 6,
    plot_type = c("enrichment", "positions"),
    assay = "peak"
) {
  plot_type <- match.arg(plot_type)
  active_network <- object@active_network
  features <- object@networks[[active_network]][[cell_type]]
  if(is.null(features)) {
    stop(sprintf("No features found for cell type: %s", cell_type))
  }

  top_peaks <- rownames(features$peaks)
  if(length(top_peaks) == 0) {
    stop(sprintf("No significant peaks found for cell type: %s", cell_type))
  }

  peaks_info <- strsplit(top_peaks, "\\.")
  peaks_coords <- sapply(peaks_info, function(x) x[2])

  peaks_parts <- strsplit(peaks_coords, "-")
  peaks_ranges <- GenomicRanges::GRanges(
    seqnames = sapply(peaks_parts, `[`, 1),
    ranges = IRanges::IRanges(
      start = as.numeric(sapply(peaks_parts, `[`, 2)),
      end = as.numeric(sapply(peaks_parts, `[`, 3))
    )
  )

  standard_peaks <- paste0(
    GenomeInfoDb::seqnames(peaks_ranges), "-",
    IRanges::start(peaks_ranges), "-",
    IRanges::end(peaks_ranges)
  )

  peak_matrix <- Seurat::GetAssayData(
    object@data,
    assay = assay,
    layer = "data"
  )
  valid_peaks <- intersect(standard_peaks, rownames(peak_matrix))

  if(length(valid_peaks) == 0) {
    stop("No valid peaks found in the assay matrix")
  }

  if(length(valid_peaks) < length(top_peaks)) {
    warning(sprintf(
      "Only %d out of %d peaks were found in the assay matrix",
      length(valid_peaks),
      length(top_peaks)
    ))
  }

  enriched.motifs <- Signac::FindMotifs(
    object = object@data,
    features = valid_peaks,
    assay = assay
  )

  if(nrow(enriched.motifs) == 0) {
    stop("No enriched motifs found")
  }

  top_motifs <- head(rownames(enriched.motifs), n_motifs)
  Signac::MotifPlot(
    object = object@data,
    motifs = top_motifs
  )
}

#' @export
FormatPeakNames <- function(peaks) {
  peaks_info <- strsplit(peaks, "\\.")
  peaks_coords <- sapply(peaks_info, function(x) x[2])
  peaks_parts <- strsplit(peaks_coords, "-")

  standard_peaks <- paste0(
    sapply(peaks_parts, `[`, 1), "-",
    sapply(peaks_parts, `[`, 2), "-",
    sapply(peaks_parts, `[`, 3)
  )

  return(standard_peaks)
}

#' @export
PlotMotifEnrichmentHeatmap <- function(
    object,
    cell_types = NULL,
    n_motifs = 10,
    min_pval = 0.05,
    assay = "peak"
) {
  all_motifs <- list()

  if(is.null(cell_types)) {
    cell_types <- names(object@csn@networks[["type_specific_features"]])
  }

  peak_matrix <- Seurat::GetAssayData(object@data, assay = assay, layer = "counts")
  active_network <- object@active_network
  for(ct in cell_types) {
    features <- object@networks[[active_network]][[ct]]
    if(is.null(features)) next

    standard_peaks <- FormatPeakNames(rownames(features$peaks))
    valid_peaks <- intersect(standard_peaks, rownames(peak_matrix))
    if(length(valid_peaks) == 0) next

    enriched.motifs <- Signac::FindMotifs(
      object = object@data,
      features = valid_peaks,
      assay = assay
    )

    if(nrow(enriched.motifs) > 0) {
      enriched.motifs$cell_type <- ct
      all_motifs[[ct]] <- enriched.motifs
    }
  }

  if(length(all_motifs) == 0) {
    stop("No enriched motifs found for any cell type")
  }

  combined_motifs <- do.call(rbind, all_motifs)
  combined_motifs$motif <- rownames(combined_motifs)

  significant_motifs <- combined_motifs[combined_motifs$pvalue < min_pval,]
  if(nrow(significant_motifs) == 0) {
    stop(sprintf("No motifs found with p-value < %g", min_pval))
  }

  top_motifs <- significant_motifs[order(significant_motifs$fold.enrichment, decreasing = TRUE),]
  top_motifs <- head(top_motifs, n_motifs)

  p <- ggplot2::ggplot(
    top_motifs,
    ggplot2::aes(
      x = stats::reorder(motif, fold.enrichment),
      y = fold.enrichment,
      fill = -log10(pvalue)
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top enriched motifs",
      x = "Motif",
      y = "Fold enrichment",
      fill = "-log10(p-value)"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c()

  return(list(
    plot = p,
    motifs = top_motifs
  ))
}
