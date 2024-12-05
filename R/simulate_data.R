#' Simulate data for testing inferCSN
#'
#' @param n_cells Number of cells to simulate
#' @param n_genes Number of genes to simulate
#' @param n_tfs Number of transcription factors
#' @param n_peaks Number of ATAC peaks to simulate
#' @param n_celltypes Number of cell types to simulate
#' @param sparsity Sparsity level of the expression matrix (0-1)
#' @param atac_sparsity Sparsity level of the ATAC matrix (0-1)
#'
#' @return A list containing simulated data
#' @export
simulate_csn_data <- function(
    n_cells = 300,
    n_genes = 100,
    n_tfs = 50,
    n_peaks = 1000,
    n_celltypes = 3,
    sparsity = 0.9,
    atac_sparsity = 0.9) {
  utils::data("motif2tf", envir = environment())

  gene_names <- unique(motif2tf$tf)
  gene_names <- gene_names[!is.na(gene_names)]

  all_genes <- sample(gene_names, n_genes)
  tfs <- all_genes[seq_len(n_tfs)]
  targets <- all_genes[(n_tfs + 1):n_genes]

  chrom_lengths <- c(
    "chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559,
    "chr4" = 190214555, "chr5" = 181538259, "chr6" = 170805979,
    "chr7" = 159345973, "chr8" = 145138636, "chr9" = 138394717,
    "chr10" = 133797422, "chr11" = 135086622, "chr12" = 133275309,
    "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189,
    "chr16" = 90338345, "chr17" = 83257441, "chr18" = 80373285,
    "chr19" = 58617616, "chr20" = 64444167, "chr21" = 46709983,
    "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415
  )

  peak_regions <- vapply(seq_len(n_peaks), function(i) {
    chr <- sample(names(chrom_lengths), 1)
    max_pos <- chrom_lengths[chr] - 1000
    peak_start <- sample.int(max_pos, 1)
    peak_end <- peak_start + sample(100:1000, 1)
    paste0(chr, "-", peak_start, "-", peak_end)
  }, character(1))

  expression_matrix <- Matrix::Matrix(
    0,
    nrow = n_cells,
    ncol = n_genes,
    sparse = TRUE
  )
  peak_matrix <- Matrix::Matrix(
    0,
    nrow = n_cells,
    ncol = n_peaks,
    sparse = TRUE
  )

  expr_base <- matrix(
    stats::rexp(n_cells * n_genes, rate = 0.5),
    nrow = n_cells
  )
  expr_noise <- matrix(
    stats::rnorm(n_cells * n_genes, mean = 0, sd = 0.1),
    nrow = n_cells
  )
  peak_base <- matrix(
    stats::rexp(n_cells * n_peaks, rate = 1),
    nrow = n_cells
  )
  peak_noise <- matrix(
    stats::rnorm(n_cells * n_peaks, mean = 0, sd = 0.05),
    nrow = n_cells
  )

  cell_names <- paste0("cell_", seq_len(n_cells))
  cell_metadata <- data.frame(
    cell = cell_names,
    celltype = paste0(
      "celltype_",
      sample(seq_len(n_celltypes), n_cells, replace = TRUE)
    ),
    row.names = cell_names,
    stringsAsFactors = FALSE
  )

  for (ct in 1:n_celltypes) {
    ct_mask <- cell_metadata$celltype == paste0("celltype_", ct)

    ct_expr <- expr_base[ct_mask, ] + expr_noise[ct_mask, ]
    ct_specific_genes <- sample(n_genes, round(n_genes * 0.3))
    ct_expr[, ct_specific_genes] <- ct_expr[, ct_specific_genes] * 5
    ct_expr[ct_expr < 0] <- 0

    ct_peak <- peak_base[ct_mask, ] + peak_noise[ct_mask, ]
    ct_specific_peaks <- sample(n_peaks, round(n_peaks * 0.4))
    ct_peak[, ct_specific_peaks] <- ct_peak[, ct_specific_peaks] * 4
    ct_peak[ct_peak < 0] <- 0

    expression_matrix[ct_mask, ] <- Matrix::Matrix(ct_expr, sparse = TRUE)
    peak_matrix[ct_mask, ] <- Matrix::Matrix(ct_peak, sparse = TRUE)
  }

  expr_idx <- sample.int(
    length(expression_matrix@x),
    round(length(expression_matrix@x) * sparsity)
  )
  peak_idx <- sample.int(
    length(peak_matrix@x),
    round(length(peak_matrix@x) * atac_sparsity)
  )

  expression_matrix@x[expr_idx] <- 0
  peak_matrix@x[peak_idx] <- 0

  colnames(expression_matrix) <- all_genes
  rownames(expression_matrix) <- cell_names
  colnames(peak_matrix) <- peak_regions
  rownames(peak_matrix) <- cell_names

  expression_matrix <- Matrix::drop0(expression_matrix)
  peak_matrix <- Matrix::drop0(peak_matrix)

  peak_annotations <- data.frame(
    peak = peak_regions,
    gene = sample(all_genes, length(peak_regions), replace = TRUE),
    distance = sample.int(10000, length(peak_regions)),
    stringsAsFactors = FALSE
  )

  return(
    list(
      expression_matrix = expression_matrix,
      peak_matrix = peak_matrix,
      tfs = tfs,
      targets = targets,
      cell_metadata = cell_metadata,
      peak_annotations = peak_annotations
    )
  )
}

#' Create a Seurat object from simulated data
#'
#' @param sim_data Output from \code{\link{simulate_csn_data}}, or a list with
#'   \code{expression_matrix}, \code{peak_matrix}, \code{tfs}, \code{targets},
#'   \code{cell_metadata}, and \code{peak_annotations}.
#' @param min_cells Include genes/peaks detected in at least this many cells
#' @param min_features Include cells where at least this many genes are detected
#' @param variable_features Number of variable features to select
#'
#' @return A Seurat object with RNA and ATAC assays
#' @export
#'
#' @examples
#' \dontrun{
#' sim_data <- simulate_csn_data()
#' seurat_obj <- create_seurat_object(sim_data)
#' }
create_seurat_object <- function(
    sim_data,
    min_cells = 3,
    min_features = 0,
    variable_features = 2000) {
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = t(sim_data$expression_matrix),
    meta.data = sim_data$cell_metadata,
    min.cells = min_cells,
    min.features = min_features,
    assay = "RNA"
  )

  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    assay = "RNA"
  )

  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = variable_features,
    assay = "RNA"
  )

  seurat_obj <- Seurat::ScaleData(
    seurat_obj,
    features = rownames(seurat_obj),
    assay = "RNA"
  )

  peak_names <- colnames(sim_data$peak_matrix)

  peak_parts <- do.call(rbind, strsplit(peak_names, "-"))
  grange_counts <- GenomicRanges::GRanges(
    seqnames = peak_parts[, 1],
    ranges = IRanges::IRanges(
      start = as.numeric(peak_parts[, 2]),
      end = as.numeric(peak_parts[, 3])
    )
  )

  standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
  grange_use <- as.vector(
    GenomeInfoDb::seqnames(grange_counts)
  ) %in% standard_chroms
  peak_matrix_filtered <- sim_data$peak_matrix[, grange_use]
  grange_counts <- grange_counts[grange_use]

  peak_names_filtered <- paste0("peak_", seq_len(sum(grange_use)))
  names(grange_counts) <- peak_names_filtered
  colnames(peak_matrix_filtered) <- peak_names_filtered

  all_genes <- unique(c(sim_data$tfs, sim_data$targets))
  gene_chrs <- paste0(
    "chr",
    sample(
      c(1:22, "X", "Y"),
      length(all_genes),
      replace = TRUE
    )
  )
  gene_starts <- sample(1:3e8, length(all_genes))
  gene_ends <- gene_starts + 1000
  gene_strands <- sample(c("+", "-"), length(all_genes), replace = TRUE)

  gene_ranges <- GenomicRanges::GRanges(
    seqnames = gene_chrs,
    ranges = IRanges::IRanges(
      start = gene_starts,
      end = gene_ends
    ),
    strand = gene_strands,
    gene_name = all_genes
  )
  GenomeInfoDb::genome(gene_ranges) <- "hg38"

  atac_counts <- t(peak_matrix_filtered)

  cell_names <- rownames(sim_data$expression_matrix)

  rownames(atac_counts) <- peak_names_filtered
  colnames(atac_counts) <- cell_names

  atac_assay <- Signac::CreateChromatinAssay(
    counts = atac_counts,
    ranges = grange_counts,
    genome = "hg38",
    fragments = NULL,
    min.cells = min_cells,
    min.features = min_features,
    annotation = gene_ranges
  )

  seurat_obj[["ATAC"]] <- atac_assay

  Seurat::DefaultAssay(seurat_obj) <- "RNA"

  return(seurat_obj)
}
