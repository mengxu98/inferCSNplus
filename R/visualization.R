##' Plot links of inferCSN on focused genes
#'
#' @param weight_table inferCSN result
#' @param gene_anno the genome annotation
#' @param marker the focused gene
#' @param cutoff the cutoff of importance scores.
#' CREs are links with importance scores higher than cutoff
#' @param upstream the number of distance upstream to TSS
#' @param downstream the number of distance downstream to TSS
#' @return ggplot2 object
#' @export
plot.connections <- function(
    weight_table,
    gene_anno = NULL,
    marker,
    cutoff = NULL,
    upstream = 250000,
    downstream = 250000) {
  if (is.null(gene_anno)) {
    stop(
      paste0("No gene annotation object, please run following code:", "
        temp <- getwd()
        if (!file.exists(paste0(temp, 'Homo_sapiens.GRCh38.100.gtf.gz'))) {
          options(timeout = 2000)
          download.file(
            'ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz',
            paste0(temp, '/Homo_sapiens.GRCh38.100.gtf.gz')
          )
        }
        gene_anno <- rtracklayer::readGFF(
          paste0(temp, '/Homo_sapiens.GRCh38.100.gtf.gz')
        )
       # rename some columns to match requirements
       gene_anno$chromosome <- paste0('chr', gene_anno$seqid)
       gene_anno$gene <- gene_anno$gene_id
       gene_anno$transcript <- gene_anno$transcript_id
       gene_anno$symbol <- gene_anno$gene_name
    ")
    )
  }
  # take out the result of marker
  if (length(which(weight_table$gene %in% marker)) > 0) {
    weight_table <- weight_table[which(weight_table$gene %in% marker), , drop = FALSE]
    conns <- weight_table[, 5:7]
    names(conns) <- c("Peak1", "Peak2", "coaccess")
    rownames(conns) <- NULL

    # normalize
    conns$coaccess <- log10(abs(as.numeric(conns$coaccess)) * 100)
    if (is.null(cutoff)) {
      cutoff <- max(
        0, max(conns$coaccess[which(weight_table$function_type == "MC")])
      )
    }
    cicero::plot_connections(
      conns,
      weight_table$Chr[1],
      as.numeric(weight_table$Starts[1]) - upstream,
      as.numeric(weight_table$Starts[1]) + downstream,
      gene_model = gene_anno,
      coaccess_cutoff = cutoff,
      collapseTranscripts = "longest",
      peak_color = "#B4656F",
      connection_color = "#1563cc",
      connection_color_legend = TRUE,
      alpha_by_coaccess = FALSE,
      connection_width = 2,
      connection_ymax = NULL,
      gene_model_color = "#116211",
      gene_model_shape = c("smallArrow", "box"),
      comparison_track = NULL,
      comparison_coaccess_cutoff = 0,
      comparison_peak_color = "#B4656F",
      comparison_connection_color = "#7F7CAF",
      comparison_connection_color_legend = TRUE,
      comparison_connection_width = 2,
      comparison_ymax = NULL,
      include_axis_track = TRUE,
      return_as_list = FALSE,
      viewpoint = NULL,
      comparison_viewpoint = TRUE,
      viewpoint_color = "#F0544F",
      viewpoint_fill = "#EFD8D7",
      viewpoint_alpha = 0.5
    )
  } else {
    message("Please try other markers!")
  }
}
