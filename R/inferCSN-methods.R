#' @include setClass.R
#' @include setGenerics.R
#' @include inferCSN.R

#' @param aggregate Logical, whether to aggregate the data.
#' @param k_neigh Number of cells to be aggregated per cluster.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary
#' @param max_overlap The maximum overlapping ratio of two clusters.
#' @param reduction_name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, should accessibility values be normalized by size factor.
#' @param genome_info Genome information.
#' @param high_corr_cutoff atacbinary
#' @param low_corr_cutoff atacbinary
#' @param rescued Logical
#'
#' @export
#' @method inferCSN Seurat
#'
#' @rdname inferCSN
#' @examples
#' \dontrun{
#' data("promoter_regions_hg38")
#' seurat_object <- inferCSN(
#'  seurat_object,
#'  genome_info = promoter_regions_hg38
#' )
#' }
setMethod(
  f = "inferCSN",
  signature = signature(object = "Seurat"),
  definition = function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    percent_samples = 1,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    cores = 1,
    verbose = FALSE,
    aggregate = TRUE,
    k_neigh = 50,
    atacbinary = TRUE,
    max_overlap = 0.8,
    reduction_name = NULL,
    size_factor_normalize = FALSE,
    genome_info = NULL,
    high_corr_cutoff = NULL,
    low_corr_cutoff = NULL,
    rescued = TRUE,
    ...) {
    if (verbose) {
      message(paste("Running start for <", class(object)[1], "object >."))
    }
    object$cluster <- Seurat::Idents(object)

    clusters <- as.character(unique(object$cluster))
    network_table_final_list <- list()
    network_table_rna_list <- list()
    network_table_atac_list <- list()
    for (c in seq_along(clusters)) {
      cluster <- clusters[c]
      if (verbose) {
        message(paste0("Running for cluster: ", cluster, "."))
      }
      object_sub <- subset(object, cluster == cluster)
      if (ncol(object_sub) < 50) {
        return()
        next
      }

      targets <- dynamic_genes(object_sub)

      if (aggregate) {
        if ("aggregated_data" %in% names(Seurat::Misc(object_sub))) {
          agg_data <- Seurat::Misc(object_sub, slot = "aggregated_data")
        } else {
          agg_data <- aggregating_data(
            object_sub,
            k_neigh = k_neigh,
            atacbinary = atacbinary,
            max_overlap = max_overlap,
            reduction_name = reduction_name,
            size_factor_normalize = size_factor_normalize,
            verbose = verbose
          )
          Seurat::Misc(object_sub, slot = "aggregated_data") <- agg_data
        }

        if ("RNA" %in% names(agg_data)) {
          data_rna <- as_matrix(agg_data$RNA)
        }

        if ("ATAC" %in% names(agg_data)) {
          data_atac <- as_matrix(agg_data$ATAC)
        }
      } else {
        if ("RNA" %in% names(object_sub@assays)) {
          data_rna <- as_matrix(object_sub@assays$RNA$data)
        }
        if ("ATAC" %in% names(object_sub@assays)) {
          data_atac <- as_matrix(object_sub@assays$ATAC$data)
        }
      }

      if ("RNA" %in% names(object_sub@assays)) {
        Seurat::DefaultAssay(object_sub) <- "RNA"
        network_table_rna <- inferCSN(
          t(data_rna),
          penalty = penalty,
          algorithm = algorithm,
          cross_validation = cross_validation,
          seed = seed,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators = regulators,
          targets = targets,
          regulators_num = regulators_num,
          verbose = verbose,
          cores = cores,
          ...
        )

        network_table_rna_list[[c]] <- network_table_rna
      }

      if ("ATAC" %in% names(object_sub@assays)) {
        rownames(data_atac) <- gsub("-", "_", rownames(data_atac))

        if (is.null(genome_info)) {
          stop("---No genome data provided.
             Please run: data('promoter_regions_hg38') or data('promoter_regions_hg19'),
             and: genome_info <- promoter_regions_hg38 or genome_info <- promoter_regions_hg19,
             then set parameter: genome_info = genome_info,
             or provide data by yourself.")
        }

        if (is.null(targets)) {
          # Obtain candidate regions of focus target genes
          if ("RNA" %in% names(object_sub@assays)) {
            Seurat::DefaultAssay(object_sub) <- "RNA"
            # targets <- Seurat::VariableFeatures(object = object_sub)
            if ("all_markers_list" %in% names(Seurat::Misc(object_sub))) {
              all_markers_list <- Seurat::Misc(object_sub, slot = "all_markers_list")
            }
            all_markers_list <- as.data.frame(all_markers_list)
            focused_markers <- all_markers_list[which(all_markers_list$cluster %in% cluster), , drop = FALSE]
            targets <- focused_markers$gene

            targets <- lapply(targets, function(x) strsplit(x, "[.]")[[1]][1])
            targets <- unique(unlist(targets))
            targets <- genome_info$genes[which(genome_info$genes %in% targets)]
            genome_info_used <- genome_info[which(genome_info$genes %in% targets), ]
          } else {
            genome_info_used <- genome_info
          }
        } else {
          targets <- lapply(targets, function(x) strsplit(x, "[.]")[[1]][1])
          targets <- unique(unlist(targets))
          targets <- genome_info$genes[which(genome_info$genes %in% targets)]
          genome_info_used <- genome_info[which(genome_info$genes %in% targets), ]
        }

        targets <- genome_info_used$genes
        chr <- genome_info_used$Chrom
        starts <- genome_info_used$Starts
        ends <- genome_info_used$Ends

        network_table_atac <- .inferCSN.atac(
          peak_matrix = data_atac,
          penalty = penalty,
          algorithm = algorithm,
          cross_validation = cross_validation,
          seed = seed,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators = targets,
          targets = targets,
          regulators_num = regulators_num,
          verbose = verbose,
          cores = cores,
          chr,
          starts,
          ends,
          high_corr_cutoff = high_corr_cutoff,
          low_corr_cutoff = low_corr_cutoff,
          rescued = rescued
        )

        network_table_atac <- network_table_atac[order(
          abs(as.numeric(network_table_atac$weight)),
          decreasing = TRUE
        ), ]
        network_table_atac_list[[c]] <- network_table_atac
        network_table_atac_sub <- extract_network(object_sub, cluster = clusters[c])
        names(network_table_atac_sub) <- c("regulator", "target", "celltype", "types")
      }

      if (all(c("RNA", "ATAC") %in% names(object_sub@assays))) {
        network_table_final <- merge(
          network_table_atac_sub,
          network_table_rna,
          by = c("regulator", "target"),
          all.x = TRUE
        )
        network_table_final <- na.omit(network_table_final)[, 1:5]
        network_table_final$weight <- network_table_final$weight / sum(abs(network_table_final$weight))
        network_table_final <- network_table_final[order(abs(as.numeric(network_table_final$weight)), decreasing = TRUE), ]
        network_table_final_list[[c]] <- network_table_final
      }
      if (verbose) message(paste0("Run done for cluster: ", cluster, "."))
    }

    if ("RNA" %in% names(object_sub@assays)) {
      names(network_table_rna_list) <- clusters
    }
    if ("ATAC" %in% names(object_sub@assays)) {
      names(network_table_atac_list) <- clusters
    }
    if (all(c("RNA", "ATAC") %in% names(object_sub@assays))) {
      names(network_table_final_list) <- clusters
    }

    # save result
    Seurat::Misc(object, slot = "network_table_rna_list") <- network_table_rna_list
    Seurat::Misc(object, slot = "network_table_atac_list") <- network_table_atac_list
    Seurat::Misc(object, slot = "network_table_final_list") <- network_table_final_list

    if (verbose) message("Run done.")

    return(object)
  }
)

.inferCSN.atac <- function(
    peak_matrix,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    percent_samples = 1,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    chr,
    starts,
    ends,
    high_corr_cutoff = NULL,
    low_corr_cutoff = NULL,
    rescued = TRUE,
    ...) {
  cores <- 1
  matrix <- as_matrix(peak_matrix)
  rm(peak_matrix)
  peaks <- rownames(matrix)

  if (!is.null(regulators)) {
    regulators_matrix <- matrix[, intersect(colnames(matrix), regulators)]
  } else {
    regulators_matrix <- matrix
  }

  if (!is.null(targets)) {
    targets_matrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targets_matrix <- matrix
  }

  cores <- min(
    (parallel::detectCores(logical = FALSE) - 1), cores, length(targets)
  )

  enhancers_matrix_list <- list()
  promoters_matrix_list <- list()
  weight_list <- list()
  if (cores == 1) {
    if (verbose) message("Using 1 core.")
    # Format progress information
    format <- "Running [:bar] :percent, No.:current of :total targets, :elapsed."
    pb <- progress::progress_bar$new(
      format = format,
      total = length(targets),
      clear = TRUE,
      width = 80
    )

    peaks <- rownames(matrix)

    for (i in seq_along(targets)) {
      if (verbose) pb$tick()

      p1 <- paste0(chr[i], ":", (starts[i] - 500), "-", starts[i])
      p2 <- paste0(chr[i], ":", (starts[i] - 250000), "-", (starts[i] + 250000))
      promoters <- find.overlapping.coordinates(peaks, p1)
      enhancers <- find.overlapping.coordinates(peaks, p2)
      enhancers <- setdiff(enhancers, promoters)

      if (length(promoters) > 0 && length(enhancers) > 1) {
        id1 <- na.omit(match(promoters, peaks))
        id2 <- na.omit(match(enhancers, peaks))
        X <- matrix[setdiff(id2, id1), ]
        enhancers_matrix_list[[i]] <- X
        Y <- matrix[id1, ]
        if (length(id1) > 1) {
          Y <- colSums(Y)
        }
        Y <- t(as_matrix(Y))
        rownames(Y) <- peaks[id1[1]]
        promoters_matrix_list[[i]] <- Y
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as_matrix(Y)
        target <- rownames(Y)
        regulators <- rownames(X)
        weight_list[[i]] <- single_network(
          matrix = t(matrix),
          regulators = regulators,
          target = target,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      }
    }
  } else {
    doParallel::registerDoParallel(cores = cores)
    if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weight_list <- foreach::foreach(
      i = seq_along(targets),
      .export = c("single_network", "sparse.regression")
    ) %dopar% {
      p1 <- paste0(chr[i], ":", (starts[i] - 500), "-", starts[i])
      p2 <- paste0(chr[i], ":", (starts[i] - 250000), "-", (starts[i] + 250000))
      promoters <- find.overlapping.coordinates(peaks, p1)
      enhancers <- find.overlapping.coordinates(peaks, p2)
      enhancers <- setdiff(enhancers, promoters)

      if (length(promoters) > 0 && length(enhancers) > 1) {
        id1 <- na.omit(match(promoters, peaks))
        id2 <- na.omit(match(enhancers, peaks))
        X <- matrix[setdiff(id2, id1), ]
        enhancers_matrix_list[[i]] <- X
        Y <- matrix[id1, ]
        if (length(id1) > 1) {
          Y <- colSums(Y)
        }
        Y <- t(as_matrix(Y))
        rownames(Y) <- peaks[id1[1]]
        promoters_matrix_list[[i]] <- Y
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as_matrix(Y)
        target <- rownames(Y)
        regulators <- rownames(X)
        single_network(
          matrix = t(matrix),
          regulators = regulators,
          target = target,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      }
    }
    doParallel::stopImplicitCluster()
  }

  network_table <- purrr::list_rbind(weight_list)
  network_table <- network_table[order(abs(as.numeric(network_table$weight)), decreasing = TRUE), ]
  # colnames(network_table) <- c("Peak1", "Peak2", "weight")

  # variable selection
  if (is.null(high_corr_cutoff)) {
    high_corr_cutoff <- max(stats::quantile(network_table$weight, 0.50), 0.001)
  }
  if (is.null(low_corr_cutoff)) {
    low_corr_cutoff <- min(0.001, stats::quantile(network_table$weight, 0.25))
  }

  for (i in 1:length(weight_list)) {
    if (!is.null(weight_list[[i]])) {
      network_table_h <- weight_list[[i]]
      coefficients <- network_table_h[, 3]
      index1 <- which(coefficients > high_corr_cutoff)
      index2 <- intersect(which(coefficients > low_corr_cutoff), which(coefficients <= high_corr_cutoff))
      index3 <- which(coefficients <= low_corr_cutoff)

      # do data frame: gene, starts, end, peak1, peak2,weight,function_type
      function_type <- rep(NA, length(coefficients))
      function_type[index1] <- "high_corr"
      function_type[index2] <- "medain_corr"
      function_type[index3] <- "low_corr"

      # rescue highly correlated CREs
      if (rescued) {
        if (i <= length(enhancers_matrix_list)) {
          X <- enhancers_matrix_list[[i]]
          Y <- promoters_matrix_list[[i]]
          CPi <- abs(cor(t(X)))
          for (p in 1:nrow(CPi)) {
            CPi[p, p] <- 0
          }
          # focus on high_corr rows
          hic_index <- which(rownames(X) %in% network_table_h[, 2][index1])
          other_index <- which(rownames(X) %in% network_table_h[, 2][-index1])
          CPi_sub <- CPi[hic_index, other_index, drop = FALSE]
          flag_matrix <- matrix(0, nrow = nrow(CPi_sub), ncol = ncol(CPi_sub))
          flag_matrix[which(CPi_sub > 0.25)] <- 1
          correlated_index <- which(colSums(flag_matrix) > 0)
          if (!is.null(correlated_index)) {
            function_type[network_table_h$Peak2 %in% rownames(X)[other_index[correlated_index]]] <- "high_corr"
          }
        }
      }
      weight_list[[i]] <- cbind(
        data.frame(
          gene = targets[i],
          Chr = chr[i],
          Starts = starts[i],
          Ends = ends[i]
        ),
        cbind(
          network_table_h,
          function_type = function_type
        )
      )
    }
  }
  network_table <- purrr::list_rbind(weight_list)
  network_table$Starts <- as.numeric(network_table$Starts)
  network_table$Ends <- as.numeric(network_table$Ends)
  network_table$weight <- abs(as.numeric(network_table$weight))

  return(network_table)
}

#' @param genes A character vector with the target genes to consider for CSN inference.
#' Takes all \code{VariableFeatures} in the object per default.
#' @param network_name network_name.
#' @param peak_to_gene_method Character specifying the method to
#' link peak overlapping motif regions to nearby genes. One of \code{Signac} or \code{GREAT}.
#' @param upstream Integer defining the distance upstream of the gene to consider as potential regulatory region.
#' @param downstream Integer defining the distance downstream of the gene to consider as potential regulatory region.
#' @param extend Integer defining the distance from the upstream and downstream of the basal regulatory region.
#' Only used of `peak_to_gene_method = 'GREAT'`.
#' @param only_tss Logical. Measure distance from the TSS (\code{TRUE}) or from the entire gene body (\code{FALSE}).
#' @param peak_to_gene_domains \code{GenomicRanges} object with regulatory regions for each gene.
#' @param tf_cor Threshold for TF - target gene correlation.
#' @param peak_cor Threshold for binding peak - target gene correlation.
#' @param aggregate_rna_col aggregate_rna_col
#' @param aggregate_peaks_col aggregate_peaks_col
#' @param method A character string indicating the method to fit the model.
#' * \code{'srm'} - Sparse Regression Model.
#' * \code{'glm'} - Generalized Liner Model with \code{\link[stats]{glm}}.
#' * \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
#' * \code{'brms'} - Bayesian Regression Models using \code{\link[brms]{brms-package}}.
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' * \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.
#' * \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[stats]{family}} for mode details.
#' @param interaction_term The interaction term to use in the model between TF and binding site.
#' * \code{'+'} for additive interaction.
#' * \code{':'} for 'multiplicative' interaction.
#' * \code{'*'} for crossing interaction, i.e. additive AND 'multiplicative'.
#' For more info, see \code{\link[stats]{formula}}
#' @param scale Logical. Whether to z-transform the expression and accessibility matrices.
#' @param adjust_method Method for adjusting p-values.
#'
#' @return A CSNObject.
#'
#' @rdname inferCSN
#' @method inferCSN CSNObject
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(
    object = "CSNObject"),
  definition = function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    percent_samples = 1,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    cores = 1,
    verbose = TRUE,
    genes = NULL,
    network_name = paste0(method, "_network"),
    peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    tf_cor = 0.1,
    peak_cor = 0.,
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c("srm", "glm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge", "bayesian_ridge"),
    alpha = 0.5,
    family = "gaussian",
    interaction_term = ":",
    adjust_method = "fdr",
    scale = FALSE,
    ...) {
    # Match args
    method <- match.arg(method)
    peak_to_gene_method <- match.arg(peak_to_gene_method)

    # Fit models
    object <- fit_models(
      object = object,
      genes = genes,
      network_name = network_name,
      peak_to_gene_method = peak_to_gene_method,
      upstream = upstream,
      downstream = downstream,
      extend = extend,
      only_tss = only_tss,
      peak_to_gene_domains = peak_to_gene_domains,
      tf_cor = tf_cor,
      peak_cor = peak_cor,
      aggregate_rna_col = aggregate_rna_col,
      aggregate_peaks_col = aggregate_peaks_col,
      method = method,
      alpha = alpha,
      family = family,
      interaction_term = interaction_term,
      adjust_method = adjust_method,
      scale = scale,
      verbose = verbose,
      cores = cores,
      ...
    )

    return(object)
  }
)

#' @inheritParams inferCSN
#'
#' @return A CSNObject.
#'
#' @rdname fit_models
#' @method fit_models CSNObject
fit_models.CSNObject <- function(
    object,
    genes = NULL,
    network_name = paste0(method, "_network"),
    peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    tf_cor = 0.1,
    peak_cor = 0.,
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c("srm", "glm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge", "bayesian_ridge"),
    interaction_term = ":",
    adjust_method = "fdr",
    scale = FALSE,
    verbose = TRUE,
    cores = 1,
    ...) {
  # Match args
  method <- match.arg(method)
  peak_to_gene_method <- match.arg(peak_to_gene_method)

  # Get variables from object
  params <- Params(object)
  motif2tf <- NetworkTFs(object)
  if (is.null(motif2tf)) {
    stop("Motif matches have not been found. Please run `find_motifs()` first.")
  }
  gene_annot <- Signac::Annotation(GetAssay(object, params$peak_assay))
  if (is.null(gene_annot)) {
    stop("Please provide a gene annotation for the ChromatinAssay.")
  }
  # Select target genes for CSN inference
  if (is.null(genes)) {
    genes <- VariableFeatures(object, assay = params$rna_assay)
    if (is.null(genes)) {
      stop("Please provide a set of features or run `FindVariableFeatures()`")
    }
  }

  # Get assay data or summary
  if (is.null(aggregate_rna_col)) {
    gene_data <- Matrix::t(LayerData(
      object,
      assay = params$rna_assay, layer = "data"
    ))
    gene_groups <- TRUE
  } else {
    gene_data <- GetAssaySummary(
      object,
      assay = params$rna_assay,
      group_name = aggregate_rna_col,
      verbose = FALSE
    )
    gene_groups <- object@meta.data[[aggregate_rna_col]]
  }

  if (is.null(aggregate_peaks_col)) {
    peak_data <- Matrix::t(LayerData(
      object,
      assay = params$peak_assay, layer = "data"
    ))
    peak_groups <- TRUE
  } else {
    peak_data <- GetAssaySummary(
      object,
      assay = params$peak_assay,
      group_name = aggregate_peaks_col,
      verbose = FALSE
    )
    peak_groups <- object@data@meta.data[[aggregate_peaks_col]]
  }

  # Select genes to use by intersecting annotated genes with all
  # detected genes in the object
  features <- intersect(gene_annot$gene_name, genes) %>%
    intersect(rownames(GetAssay(object, params$rna_assay)))
  gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]

  # Get regions
  regions <- NetworkRegions(object)
  peak_data <- peak_data[, regions@peaks]
  colnames(peak_data) <- rownames(regions@motifs@data)
  peaks2motif <- regions@motifs@data

  # Find candidate regions near gene bodies
  if (is.null(peak_to_gene_domains)) {
    log_message("Selecting candidate regulatory regions near genes", verbose = verbose)
    peaks_near_gene <- find_peaks_near_genes(
      peaks = regions@ranges,
      method = peak_to_gene_method,
      genes = gene_annot,
      upstream = upstream,
      downstream = downstream,
      only_tss = only_tss
    )
  } else {
    log_message("Selecting candidate regulatory regions in provided domains", verbose = verbose)
    peaks_near_gene <- find_peaks_near_genes(
      peaks = regions@ranges,
      method = "Signac",
      genes = peak_to_gene_domains,
      upstream = 0,
      downstream = 0,
      only_tss = FALSE
    )
  }

  peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups = colnames(peaks_near_gene), fun = "sum")

  # Select peaks passing criteria
  peaks_at_gene <- as.logical(colMaxs(peaks2gene))
  peaks_with_motif <- as.logical(rowMaxs(peaks2motif * 1))

  # Subset data to good peaks
  peaks_use <- peaks_at_gene & peaks_with_motif
  peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
  peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
  peak_data <- peak_data[, peaks_use, drop = FALSE]

  log_message("Preparing model input", verbose = verbose)
  tfs_use <- colnames(motif2tf)
  motif2tf <- motif2tf[, tfs_use, drop = FALSE]

  log_message("Fitting models for ", length(features), " target genes", verbose = verbose)
  # Loop through features and fit models/run CV for each
  names(features) <- features
  model_fits <- parallelize_fun(
    features, function(g) {
      # Select peaks near gene
      if (!g %in% rownames(peaks2gene)) {
        log_message("Warning: ", g, " not found in EnsDb", verbose = verbose == 2)
        return()
      }
      gene_peaks <- as.logical(peaks2gene[g, ])
      if (sum(gene_peaks) == 0) {
        log_message("Warning: No peaks found near ", g, verbose = verbose == 2)
        return()
      }

      # Select peaks correlating with target gene expression
      g_x <- gene_data[gene_groups, g, drop = FALSE]
      peak_x <- peak_data[peak_groups, gene_peaks, drop = FALSE]
      peak_g_cor <- as(sparse_cor(peak_x, g_x), "generalMatrix")
      peak_g_cor[is.na(peak_g_cor)] <- 0
      peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > peak_cor]
      if (length(peaks_use) == 0) {
        log_message("Warning: No correlating peaks found for ", g, verbose = verbose == 2)
        return()
      }
      peak_x <- peak_x[, peaks_use, drop = FALSE]
      peak_motifs <- peaks2motif[gene_peaks, , drop = FALSE][peaks_use, , drop = FALSE]

      # Select TFs with motifs in peaks
      gene_peak_tfs <- purrr::map(rownames(peak_motifs), function(p) {
        x <- as.logical(peak_motifs[p, ])
        peak_tfs <- colMaxs(motif2tf[x, , drop = FALSE])
        peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
        peak_tfs <- setdiff(peak_tfs, g)
        return(peak_tfs)
      })
      names(gene_peak_tfs) <- rownames(peak_motifs)

      # Check correlation of peaks with target gene
      gene_tfs <- purrr::reduce(gene_peak_tfs, union)
      tf_x <- gene_data[gene_groups, gene_tfs, drop = FALSE]
      tf_g_cor <- as(sparse_cor(tf_x, g_x), "generalMatrix")
      tf_g_cor[is.na(tf_g_cor)] <- 0
      tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > tf_cor]
      if (length(tfs_use) == 0) {
        log_message("Warning: No correlating TFs found for ", g, verbose = verbose == 2)
        return()
      }
      tf_g_corr_df <- tibble::as_tibble(
        tf_g_cor[unique(tfs_use), , drop = FALSE],
        rownames = "tf",
        .name_repair = "check_unique"
      ) %>%
        rename("tf" = 1, "corr" = 2)

      # Filter TFs and make formula string
      frml_string <- purrr::map(names(gene_peak_tfs), function(p) {
        peak_tfs <- gene_peak_tfs[[p]]
        peak_tfs <- peak_tfs[peak_tfs %in% tfs_use]
        if (length(peak_tfs) == 0) {
          return()
        }
        peak_name <- stringr::str_replace_all(p, "-", "_")
        tf_name <- stringr::str_replace_all(peak_tfs, "-", "_")
        formula_str <- paste(
          paste(peak_name, interaction_term, tf_name, sep = " "),
          collapse = " + "
        )
        return(list(tfs = peak_tfs, frml = formula_str))
      })
      frml_string <- frml_string[!purrr::map_lgl(frml_string, is.null)]
      if (length(frml_string) == 0) {
        log_message("Warning: No valid peak:TF pairs found for ", g, verbose = verbose == 2)
        return()
      }

      target <- stringr::str_replace_all(g, "-", "_")
      model_frml <- stats::as.formula(
        paste0(target, " ~ ", paste0(purrr::map(frml_string, function(x) x$frml), collapse = " + "))
      )

      # Get expression data
      nfeats <- sum(purrr::map_dbl(frml_string, function(x) length(x$tfs)))
      gene_tfs <- purrr::reduce(purrr::map(frml_string, function(x) x$tfs), union)
      gene_x <- gene_data[gene_groups, union(g, gene_tfs), drop = FALSE]
      model_mat <- as.data.frame(cbind(gene_x, peak_x))
      if (scale) model_mat <- as.data.frame(scale(as_matrix(model_mat)))
      colnames(model_mat) <- stringr::str_replace_all(colnames(model_mat), "-", "_")

      log_message("Fitting model with ", nfeats, " variables for ", g, verbose = verbose == 2)
      result <- try(fit_model(
        model_frml,
        data = model_mat,
        method = method,
        ...
      ), silent = TRUE)
      if (any(class(result) == "try-error")) {
        log_message("Warning: Fitting model failed for ", g, verbose = verbose)
        log_message(result, verbose = verbose == 2)
        return()
      } else {
        result$gof$nvariables <- nfeats
        result$corr <- tf_g_corr_df
        return(result)
      }
    }, verbose = verbose, cores = cores)

  model_fits <- model_fits[!purrr::map_lgl(model_fits, is.null)]
  if (length(model_fits) == 0) {
    log_message("Warning: Fitting model failed for all genes.", verbose = verbose)
  }

  coefs <- purrr::map_dfr(model_fits, function(x) x$coefs, .id = "target")
  coefs <- format_coefs(
    coefs,
    term = interaction_term,
    adjust_method = adjust_method
  )
  corrs <- purrr::map_dfr(model_fits, function(x) x$corr, .id = "target")
  if (nrow(coefs) > 0) {
    coefs <- suppressMessages(left_join(coefs, corrs))
  }
  gof <- purrr::map_dfr(model_fits, function(x) x$gof, .id = "target")

  params <- list()
  params[["method"]] <- method
  params[["family"]] <- family
  params[["dist"]] <- c("upstream" = upstream, "downstream" = downstream)
  params[["only_tss"]] <- only_tss
  params[["interaction"]] <- interaction_term
  params[["tf_cor"]] <- tf_cor
  params[["peak_cor"]] <- peak_cor

  network_obj <- methods::new(
    Class = "Network",
    features = features,
    coefs = coefs,
    fit = gof,
    params = params
  )
  object@csn@networks[[network_name]] <- network_obj
  object@csn@active_network <- network_name

  return(object)
}

#' @param clusters clusters.
#' @return A CSNObject_list.
#'
#' @rdname inferCSN
#' @method inferCSN CSNObjectList
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(object = "CSNObjectList"),
  definition = function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    percent_samples = 1,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    cores = 1,
    verbose = TRUE,
    genes = NULL,
    network_name = paste0(method, "_network"),
    peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    tf_cor = 0.1,
    peak_cor = 0.,
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c("srm", "glm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge", "bayesian_ridge"),
    alpha = 0.5,
    family = "gaussian",
    interaction_term = ":",
    adjust_method = "fdr",
    scale = FALSE,
    clusters = NULL,
    ...) {
    method <- match.arg(method)
    peak_to_gene_method <- match.arg(peak_to_gene_method)

    if (is.null(clusters)) {
      clusters <- names(object@data)
    }

    objects_list <- purrr::map(
      clusters, function(x) {
        log_message("Running for cluster: ", x, verbose = verbose)

        inferCSN(
          object = object@data[[x]],
          genes = genes,
          network_name = network_name,
          peak_to_gene_method = peak_to_gene_method,
          upstream = upstream,
          downstream = downstream,
          extend = extend,
          only_tss = only_tss,
          peak_to_gene_domains = peak_to_gene_domains,
          tf_cor = tf_cor,
          peak_cor = peak_cor,
          aggregate_rna_col = aggregate_rna_col,
          aggregate_peaks_col = aggregate_peaks_col,
          method = method,
          alpha = alpha,
          family = family,
          interaction_term = interaction_term,
          adjust_method = adjust_method,
          scale = scale,
          verbose = verbose,
          cores = cores,
          ...
        )
      }
    )
    names(objects_list) <- clusters

    return(objects_list)
  }
)
