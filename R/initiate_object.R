#' @param object A Seurat object containing gene expression and/or chromatin accessibility data.
#' @param regulators Regulators to consider for CSN inference.
#' @param targets Targets to consider for CSN inference.
#'
#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "matrix",
  definition = function(object,
                        regulators = character(0),
                        targets = character(0),
                        ...) {
    if (length(regulators) == 0) {
      regulators <- colnames(object)
    } else {
      regulators <- intersect(
        regulators, colnames(object)
      )
    }
    if (length(targets) == 0) {
      targets <- colnames(object)
    } else {
      targets <- intersect(
        targets, colnames(object)
      )
    }

    object <- methods::new(
      Class = "Network",
      data = object,
      regulators = regulators,
      targets = targets
    )

    return(object)
  }
)

#' @param object A Seurat object containing gene expression and/or chromatin accessibility data.
#' @param celltype Selected celltype(s) to analyze. If NULL, all cells will be used.
#' @param celltype_by Column name in Seurat meta.data for cell type information.
#' If NULL, all cells will be treated as one group.
#' @param filter_mode Filter mode for identifying cell-type specific features.
#' @param filter_by Character. When filter_mode is "variable",
#' specify whether to calculate variable features using "aggregate" cells or by "celltype".
#' Default is "aggregate".
#' @param rna_assay A character vector indicating the name of the gene expression
#' assay in the \code{Seurat} object.
#' @param rna_min_pct Minimum percentage of cells expressing the gene in either population.
#' @param rna_logfc_threshold Minimum log fold-change threshold for RNA data.
#' @param rna_test_method Statistical test to use for gene markers identification.
#' @param n_variable_genes Number of variable features to select when not comparing between groups.
#' @param peak_assay A character vector indicating the name of the chromatin
#' accessibility assay in the \code{Seurat} object.
#' If NULL, only RNA-seq data will be processed.
#' @param peak_min_pct Minimum percentage of cells expressing the peak in either population.
#' @param peak_logfc_threshold Minimum log fold-change threshold for peak data.
#' @param peak_test_method Statistical test to use for peak markers identification.
#' @param n_variable_peaks Number of variable peaks to select ("q5" by default).
#' @param regions Candidate regions to consider for binding site inference.
#' If \code{NULL}, all peaks regions are considered.
#' \code{regions} could be a \code{GRanges} object,
#' or a data frame with three columns: chrom, start, end.
#' @param exclude_exons Logical. Whether to consider exons for binding site inference.
#' @param only_pos Only return positive markers.
#' @param verbose Print progress messages.
#' @param p_value Significance level threshold for filtering features (default: 0.05).
#' @param ... Additional arguments passed to marker detection functions.
#'
#' @return CSNObject object.
#'
#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "Seurat",
  definition = function(object,
                        celltype = NULL,
                        celltype_by = NULL,
                        filter_mode = c("variable", "celltype", "unfiltered"),
                        filter_by = c("aggregate", "celltype"),
                        rna_assay = "RNA",
                        rna_min_pct = 0.1,
                        rna_logfc_threshold = 0.25,
                        rna_test_method = "wilcox",
                        n_variable_genes = 2000,
                        peak_assay = NULL,
                        peak_min_pct = 0.05,
                        peak_logfc_threshold = 0.1,
                        peak_test_method = "LR",
                        n_variable_peaks = "q5",
                        regions = NULL,
                        exclude_exons = TRUE,
                        only_pos = TRUE,
                        verbose = TRUE,
                        p_value = 0.05,
                        ...) {
    if (!is.null(celltype_by) && celltype_by != "aggregate") {
      if (!celltype_by %in% colnames(object@meta.data)) {
        stop(sprintf("Column '%s' not found in object metadata", celltype_by))
      }
      Seurat::Idents(object) <- celltype_by
    } else {
      Seurat::Idents(object) <- "aggregate"
    }

    if (!is.null(celltype)) {
      if (!all(celltype %in% Seurat::Idents(object))) {
        stop("Some specified clusters not found in data")
      }
      object <- subset(object, idents = celltype)
      if (length(celltype) == 1) {
        Seurat::Idents(object) <- celltype
      }
    }

    filter_mode <- match.arg(filter_mode)
    filter_by <- match.arg(filter_by)
    celltypes <- unique(Seurat::Idents(object))

    Seurat::DefaultAssay(object) <- rna_assay
    genes_markers <- process_features(
      object = object,
      celltypes = celltypes,
      mode = filter_mode,
      by = filter_by,
      is_peak = FALSE,
      params = list(
        assay = rna_assay,
        min_pct = rna_min_pct,
        logfc_threshold = rna_logfc_threshold,
        test_method = rna_test_method,
        n_features = n_variable_genes,
        only_pos = only_pos,
        p_value = p_value
      ),
      verbose = verbose,
      ...
    )

    if (is.null(peak_assay)) {
      log_message(
        "No peak assay found.
        Please specify the peak assay with 'peak_assay' in 'initiate_object'.",
        verbose = verbose,
        message_type = "warning"
      )
    }
    peaks_markers <- NULL
    if (!is.null(peak_assay)) {
      Seurat::DefaultAssay(object) <- peak_assay
      object$nCount_peaks <- colSums(
        Seurat::GetAssayData(object, assay = peak_assay, layer = "data") > 0
      )

      peaks_markers <- process_features(
        object = object,
        celltypes = celltypes,
        mode = filter_mode,
        by = filter_by,
        is_peak = TRUE,
        params = list(
          assay = peak_assay,
          min_pct = peak_min_pct,
          logfc_threshold = peak_logfc_threshold,
          test_method = peak_test_method,
          n_features = n_variable_peaks,
          only_pos = only_pos,
          latent.vars = "nCount_peaks",
          p_value = p_value
        ),
        verbose = verbose,
        ...
      )
    }

    attributes <- process_attributes(
      object = object,
      celltypes = celltypes,
      genes_markers = genes_markers,
      peaks_markers = peaks_markers,
      verbose = verbose,
      p_value = p_value
    )

    regions_obj <- process_regions(
      object = object,
      regions = regions,
      peak_assay = peak_assay,
      exclude_exons = exclude_exons,
      verbose = verbose,
      ...
    )

    params <- list(
      peak_assay = peak_assay,
      rna_assay = rna_assay,
      exclude_exons = exclude_exons,
      filter_mode = filter_mode,
      filter_by = filter_by
    )

    print_summary(
      attributes,
      celltypes,
      peak_assay,
      verbose
    )

    object <- methods::new(
      Class = "CSNObject",
      data = object,
      metadata = list(
        celltypes = celltypes,
        n_celltypes = length(celltypes),
        attributes = attributes,
        summary = create_summary(
          attributes,
          celltypes,
          peak_assay
        )
      ),
      regions = regions_obj,
      params = params
    )

    return(object)
  }
)

#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "CSNObject",
  definition = function(object, ...) {
    initiate_object(object@data, ...)
  }
)

process_features <- function(
    object,
    celltypes,
    mode,
    by,
    is_peak,
    params,
    verbose,
    ...) {
  switch(
    EXPR = mode,
    "celltype" = process_celltype_features(
      object,
      celltypes,
      is_peak,
      params,
      verbose,
      ...
    ),
    "variable" = process_variable_features(
      object,
      celltypes,
      by,
      is_peak,
      params,
      verbose,
      ...
    ),
    "unfiltered" = process_all_features(
      object,
      celltypes,
      is_peak,
      params
    )
  )
}

process_celltype_features <- function(
    object,
    celltypes,
    is_peak,
    params,
    verbose,
    ...) {
  if (length(celltypes) <= 1) {
    return(
      process_variable_features(
        object, celltypes, "aggregate", is_peak, params, verbose, ...
      )
    )
  }

  feature_type <- if (is_peak) "peaks" else "genes"
  log_message(
    "Finding celltype-specific ", feature_type,
    verbose = verbose
  )

  purrr::map_dfr(
    celltypes, function(x) {
      markers <- suppressWarnings(
        Seurat::FindMarkers(
          object = object,
          ident.1 = x,
          assay = params$assay,
          min.pct = params$min_pct,
          logfc.threshold = params$logfc_threshold,
          test.use = params$test_method,
          only.pos = params$only_pos,
          verbose = FALSE,
          ...
        )
      )

      if (!is.null(markers) && nrow(markers) > 0) {
        markers$celltype <- x
        markers[[if (is_peak) "peak" else "gene"]] <- rownames(markers)
      }

      return(markers)
    }
  )
}

process_variable_features <- function(
    object,
    celltypes,
    by,
    is_peak,
    params,
    verbose,
    ...) {
  if (by == "aggregate") {
    object <- find_variable_features(
      object,
      is_peak,
      params
    )
    var_features <- Seurat::VariableFeatures(object)
    res <- create_feature_matrix(
      var_features,
      celltypes,
      is_peak = is_peak
    )
    return(res)
  }

  res <- purrr::map_dfr(
    celltypes, function(x) {
      cell_subset <- subset(object, idents = x)
      cell_subset <- find_variable_features(
        cell_subset,
        is_peak,
        params
      )
      var_features <- Seurat::VariableFeatures(cell_subset)
      create_feature_matrix(var_features, x, is_peak = is_peak)
    }
  )
  return(res)
}

process_all_features <- function(
    object,
    celltypes,
    is_peak,
    params,
    ...) {
  features <- rownames(
    Seurat::GetAssayData(
      object,
      assay = params$assay,
      layer = "data"
    )
  )
  res <- create_feature_matrix(features, celltypes, is_peak = is_peak)
  return(res)
}

find_variable_features <- function(
    object,
    is_peak,
    params) {
  if (is_peak) {
    Signac::FindTopFeatures(
      object,
      min.cutoff = params$n_features,
      verbose = FALSE
    )
  } else {
    Seurat::FindVariableFeatures(
      object,
      nfeatures = params$n_features,
      verbose = FALSE
    )
  }
}

.create_feature_df <- function(
    features,
    celltype,
    infinite_logfc = TRUE,
    is_peak = FALSE) {
  df <- data.frame(
    avg_log2FC = rep(if (infinite_logfc) Inf else 1, length(features)),
    p_val_adj = rep(0, length(features)),
    celltype = celltype,
    stringsAsFactors = FALSE
  )

  if (is_peak) {
    df$peak <- features
  } else {
    df$gene <- features
  }

  return(df)
}

create_feature_matrix <- function(
    features,
    celltypes,
    is_peak = FALSE) {
  if (length(celltypes) == 1) {
    celltypes <- list(celltypes)
  }
  purrr::map_dfr(
    celltypes,
    ~ .create_feature_df(features, .x, is_peak = is_peak)
  )
}

process_attributes <- function(
    object,
    celltypes,
    genes_markers,
    peaks_markers,
    verbose = TRUE,
    p_value = 0.05) {
  res <- purrr::map(
    celltypes, function(x) {
      log_message(
        "Processing results for ", x,
        verbose = verbose
      )

      sig_genes <- filter_significant_features(
        genes_markers, x, p_value
      )
      cells <- colnames(object)[Seurat::Idents(object) == x]
      sig_peaks <- if (!is.null(peaks_markers)) {
        filter_significant_features(peaks_markers, x, p_value)
      } else {
        NULL
      }

      list(
        celltype = x,
        cells = cells,
        n_cells = length(cells),
        genes = sig_genes,
        n_genes = if (is.null(sig_genes)) 0 else nrow(sig_genes),
        peaks = sig_peaks,
        n_peaks = if (is.null(sig_peaks)) 0 else nrow(sig_peaks)
      )
    }
  ) |>
    purrr::set_names(celltypes)

  return(res)
}

filter_significant_features <- function(
    markers,
    celltype,
    p_value = 0.05) {
  sig_features <- markers[
    markers$celltype == celltype &
      markers$p_val_adj < p_value,
  ]
  if (!is.null(sig_features) && nrow(sig_features) > 0) {
    sig_features[order(sig_features$avg_log2FC, decreasing = TRUE), ]
  } else {
    sig_features
  }

  return(sig_features)
}

process_regions <- function(
    object,
    regions = NULL,
    peak_assay = NULL,
    exclude_exons = TRUE,
    verbose = TRUE,
    ...) {
  if (is.null(peak_assay)) {
    res <- methods::new(
      Class = "Regions",
      ranges = GenomicRanges::GRanges(
        seqnames = character(0),
        ranges = IRanges::IRanges(
          start = integer(0),
          end = integer(0)
        )
      ),
      peaks = numeric(0),
      motifs = NULL
    )

    return(res)
  }

  log_message(
    "Processing candidate regions",
    verbose = verbose
  )

  gene_annot <- Signac::Annotation(object[[peak_assay]])
  if (is.null(gene_annot)) {
    stop("Please provide a gene annotation for the ChromatinAssay.")
  }

  peak_ranges <- Signac::StringToGRanges(
    rownames(Seurat::GetAssay(object, assay = peak_assay))
  )

  if (!is.null(regions)) {
    if (is.data.frame(regions)) {
      regions <- GenomicRanges::GRanges(
        seqnames = regions[, 1],
        ranges = IRanges::IRanges(
          start = regions[, 2],
          end = regions[, 3]
        )
      )
    }
    if (!inherits(regions, "GRanges")) {
      log_message(
        "Regions must be a GRanges object or a data frame with three columns: chrom, start, end.",
        message_type = "error"
      )
    }

    cand_olaps <- IRanges::findOverlaps(regions, peak_ranges)
    cand_ranges <- IRanges::pintersect(
      peak_ranges[S4Vectors::subjectHits(cand_olaps)],
      regions[S4Vectors::queryHits(cand_olaps)]
    )
  } else {
    cand_ranges <- peak_ranges
  }

  if (exclude_exons) {
    exon_ranges <- gene_annot[gene_annot$type == "exon", ]
    names(exon_ranges@ranges) <- NULL
    exon_ranges <- IRanges::intersect(
      exon_ranges,
      exon_ranges
    )
    exon_ranges <- GenomicRanges::GRanges(
      seqnames = exon_ranges@seqnames,
      ranges = exon_ranges@ranges
    )
    cand_ranges <- GenomicRanges::subtract(
      cand_ranges,
      exon_ranges,
      ignore.strand = TRUE
    ) |> unlist()
  }

  peak_overlaps <- IRanges::findOverlaps(
    cand_ranges,
    peak_ranges
  )
  peak_matches <- S4Vectors::subjectHits(peak_overlaps)

  res <- methods::new(
    Class = "Regions",
    ranges = cand_ranges,
    peaks = peak_matches,
    motifs = NULL
  )

  return(res)
}

print_summary <- function(
    attributes,
    celltypes,
    peak_assay,
    verbose) {
  log_message(
    "Summary of cell-type specific features:",
    verbose = verbose
  )
  for (ct in celltypes) {
    attr <- attributes[[ct]]
    if (!is.null(peak_assay)) {
      log_message(
        sprintf(
          "   %s: %d cells, %d genes, %d peaks",
          as.character(ct), attr$n_cells, attr$n_genes, attr$n_peaks
        ),
        verbose = verbose,
        cli_model = FALSE,
        timestamp = FALSE
      )
    } else {
      log_message(
        sprintf(
          "   %s: %d cells, %d genes",
          as.character(ct), attr$n_cells, attr$n_genes
        ),
        verbose = verbose,
        cli_model = FALSE,
        timestamp = FALSE
      )
    }
  }
}

create_summary <- function(
    attributes,
    celltypes,
    peak_assay) {
  res <- do.call(
    rbind,
    lapply(
      celltypes,
      function(x) {
        attr <- attributes[[x]]
        if (!is.null(peak_assay)) {
          data.frame(
            celltype = x,
            cells = attr$n_cells,
            genes = attr$n_genes,
            peaks = attr$n_peaks
          )
        } else {
          data.frame(
            celltype = x,
            cells = attr$n_cells,
            genes = attr$n_genes
          )
        }
      }
    )
  )

  return(res)
}
