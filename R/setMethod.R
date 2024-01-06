#' @param object An expression matrix, cells by genes
#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0" and "L0L2".
#' For high-dimensional and sparse data, such as single-cell sequencing data, "L0L2" is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param cross_validation Check whether cross validation is used.
#' @param n_folds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param k_folds The number of folds for sample split.
#' @param r_threshold r_threshold.
#' @param regulators Regulator genes.
#' @param targets Target genes.
#' @param regulators_num The number of non-zore coef, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param verbose Print detailed information.
#' @param cores CPU cores.
#'
#' @export
#' @method inferCSN default
#'
#' @rdname inferCSN
#'
inferCSN.default <- function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
  if (verbose) message(paste("Running start for <", class(object)[1], ">."))
  matrix <- object; rm(object)
  # Check input parameters
  check.parameters(
    matrix = matrix,
    penalty = penalty,
    algorithm = algorithm,
    cross_validation = cross_validation,
    seed = seed,
    n_folds = n_folds,
    k_folds = k_folds,
    r_threshold = r_threshold,
    regulators = regulators,
    targets = targets,
    regulators_num = regulators_num,
    verbose = verbose,
    cores = cores
  )

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
  targets <- colnames(targets_matrix)
  rm(matrix)

  cores <- min((parallel::detectCores(logical = FALSE) - 1), cores, length(targets))
  if (cores == 1) {
    if (verbose) message("Using 1 core.")
    # Format progress information
    format <- "Running [:bar] :percent, No.:current of :total genes, :elapsed."
    pb <- progress::progress_bar$new(
      format = format,
      total = length(targets),
      clear = TRUE,
      width = 80
    )

    weight_table <- purrr::map_dfr(targets, function(x) {
      if (verbose) pb$tick()
      sub.inferCSN(
        regulators_matrix = regulators_matrix,
        targets_matrix = targets_matrix,
        target = x,
        cross_validation = cross_validation,
        seed = seed,
        penalty = penalty,
        algorithm = algorithm,
        n_folds = n_folds,
        k_folds = k_folds,
        r_threshold = r_threshold,
        regulators_num = regulators_num,
        verbose = verbose
      )
    })
  } else {
    doParallel::registerDoParallel(cores = cores)
    if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")
    target <- NULL
    "%dopar%" <- foreach::"%dopar%"
    weight_list <- foreach::foreach(
      target = targets,
      .export = c("sub.inferCSN", "sparse.regression")
    ) %dopar% {
      sub.inferCSN(
        regulators_matrix = regulators_matrix,
        targets_matrix = targets_matrix,
        target = target,
        cross_validation = cross_validation,
        seed = seed,
        penalty = penalty,
        algorithm = algorithm,
        n_folds = n_folds,
        k_folds = k_folds,
        r_threshold = r_threshold,
        regulators_num = regulators_num,
        verbose = verbose
      )
    }
    weight_table <- purrr::list_rbind(weight_list)
    doParallel::stopImplicitCluster()
  }

  weight_table <- weight_table[order(abs(as.numeric(weight_table$weight)), decreasing = TRUE), ]
  if (verbose) message("Run done.")
  return(weight_table)
}

#' @export
#' @method inferCSN data.frame
#'
#' @rdname inferCSN
#'
inferCSN.data.frame <- function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
  if (verbose) message(paste("Running start for <", class(object)[1], ">."))
  if (verbose) message("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(object)
  inferCSN(
    matrix,
    penalty = penalty,
    algorithm = algorithm,
    cross_validation = cross_validation,
    seed = seed,
    n_folds = n_folds,
    k_folds = k_folds,
    r_threshold = r_threshold,
    regulators = regulators,
    targets = targets,
    regulators_num = regulators_num,
    verbose = verbose,
    cores = cores,
    ...
  )
}

#' @param aggregate aggregate
#' @param peakcalling peakcalling
#' @param macs2.path macs2.path
#' @param fragments fragments
#' @param k_neigh Number of cells to be aggregated per cluster.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary
#' @param max_overlap The maximum overlapping ratio of two clusters.
#' @param reduction_name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, should accessibility values be normalized by size factor
#' @param genome_info atacbinary
#' @param high_corr_cutoff atacbinary
#' @param low_corr_cutoff atacbinary
#' @param rescued atacbinary
#' @param ... atacbinary
#'
#' @export
#' @method inferCSN Seurat
#' @import Seurat
#' @importFrom stats cor na.omit
#' @importFrom utils data
#'
#' @rdname inferCSN
inferCSN.Seurat <- function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    aggregate = TRUE,
    peakcalling = FALSE,
    macs2.path = NULL,
    fragments = NULL,
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
  if (verbose) message(paste("Running start for <", class(object)[1], "object >."))
  object$cluster <- Seurat::Idents(object)

  # # step 0. Peak calling
  # if (("ATAC" %in% names(object@assays)) && peakcalling) {
  #   if (verbose) {
  #     message("Calling Peak")
  #   }
  #   object$cluster <- Seurat::Idents(object)
  #   if (is.null(macs2.path)) {
  #     message("Please give the path to macs2!") # https://macs3-project.github.io/MACS/
  #   }
  #   # DefaultAssay(combined) <- "ATAC"
  #   peaks <- Signac::CallPeaks(
  #     object = object,
  #     group.by = "cluster"
  #   )
  #
  #   if (is.null(fragments)) {
  #     message("Please input fragments!")
  #   }
  #   new_atac_data <- Signac:::FeatureMatrix(
  #     fragments = fragments, # fragments of original, we need to do aggregation
  #     features = peaks
  #   )
  #   object@assays$ATAC@counts <- new_atac_data
  #   if (verbose) {
  #     message("Peak calling finished")
  #   }
  # }

  object_raw <- object
  clusters <- as.character(unique(object$cluster))
  weight_table_final_list <- list()
  for (c in seq_along(clusters)) {
    cluster <- clusters[c]
    if (verbose) message(paste0("Running for cluster: ", cluster, "."))
    object <- subset(object_raw, cluster == cluster[1])

    if (aggregate) {
      if ("aggregated_data" %in% names(Seurat::Misc(object))) {
        agg_data <- Seurat::Misc(object, slot = "aggregated_data")
      } else {
        agg_data <- aggregate(
          object,
          k_neigh = k_neigh,
          atacbinary = atacbinary,
          max_overlap = max_overlap,
          reduction_name = NULL,
          size_factor_normalize = size_factor_normalize
        )
        Seurat::Misc(object, slot = "aggregated_data") <- agg_data
      }

      if ("RNA" %in% names(agg_data)) {
        data_rna <- as.matrix(agg_data$RNA)
      }

      if ("ATAC" %in% names(agg_data)) {
        data_atac <- as.matrix(agg_data$ATAC)
      }
    } else {
      if ("RNA" %in% names(object@assays)) {
        data_rna <- matrix(
          0,
          nrow = nrow(object@assays$RNA$counts),
          ncol = 1
        )
      }
      if ("ATAC" %in% names(object@assays)) {
        data_atac <- matrix(
          0,
          nrow = nrow(object@assays$ATAC@counts),
          ncol = 1
        )
      }
    }

    if ("RNA" %in% names(agg_data)) {
      Seurat::DefaultAssay(object) <- "RNA"
      if (is.null(targets)) {
        # targets <- Seurat::VariableFeatures(object = object)
        if ("all_markers_list" %in% names(Seurat::Misc(object))) {
          all_markers_list <- Seurat::Misc(object, slot = "all_markers_list")
        }
        all_markers_list <- as.data.frame(all_markers_list)
        focused_markers <- all_markers_list[which(all_markers_list$cluster %in% cluster), , drop = FALSE]
        targets <- focused_markers$gene
      }
      weight_table_rna <- inferCSN(
        t(data_rna),
        penalty = penalty,
        algorithm = algorithm,
        cross_validation = cross_validation,
        seed = seed,
        n_folds = n_folds,
        k_folds = k_folds,
        r_threshold = r_threshold,
        regulators = targets,
        targets = targets,
        regulators_num = regulators_num,
        verbose = verbose,
        cores = cores
      )

      weight_table_rna <- weight_table_rna[order(
        abs(as.numeric(weight_table_rna$weight)), decreasing = TRUE
      ), ]
      Misc(object, slot = "weight_table_rna") <- weight_table_rna
    }

    # data_atac <- peaks.filter(data_atac)
    if ("ATAC" %in% names(object@assays)) {
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
        if ("RNA" %in% names(object@assays)) {
          Seurat::DefaultAssay(object) <- "RNA"
          # targets <- Seurat::VariableFeatures(object = object)
          if ("all_markers_list" %in% names(Seurat::Misc(object))) {
            all_markers_list <- Seurat::Misc(object, slot = "all_markers_list")
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

      weight_table_atac <- .inferCSN.atac(
        peak_matrix = data_atac,
        penalty = penalty,
        algorithm = algorithm,
        cross_validation = cross_validation,
        seed = seed,
        n_folds = n_folds,
        k_folds = k_folds,
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

      weight_table_atac <- weight_table_atac[order(
        abs(as.numeric(weight_table_atac$weight)), decreasing = TRUE
      ), ]
      Seurat::Misc(object, slot = "weight_table_atac") <- weight_table_atac
      weight_table_atac_sub <- final.network(object, cluster = clusters[c])
      names(weight_table_atac_sub) <- c("regulator", "target", "celltype", "types")
    }

    if ("RNA" %in% names(object@assays)) {
      weight_table_final <- merge(
        weight_table_atac_sub,
        weight_table_rna,
        by = c("regulator", "target"),
        all.x = TRUE
      )
      weight_table_final <- na.omit(weight_table_final)[, 1:5]
      weight_table_final$weight <- weight_table_final$weight / sum(abs(weight_table_final$weight))
      weight_table_final <- weight_table_final[order(abs(as.numeric(weight_table_final$weight)), decreasing = TRUE), ]
    }
    weight_table_final_list[[c]] <- weight_table_final
    if (verbose) message(paste0("Run done for cluster: ", cluster, "."))
  }

  names(weight_table_final_list) <- clusters

  # save result
  Seurat::Misc(object_raw, slot = "weight_table") <- weight_table_final_list

  if (verbose) message("Run done.")

  return(object_raw)
}

.inferCSN.atac <- function(
    peak_matrix,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
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
  matrix <- as.matrix(peak_matrix)
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
    format <- "Running [:bar] :percent, No.:current of :total peaks, :elapsed."
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
      promoters <- cicero::find_overlapping_coordinates(peaks, p1)
      enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
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
        Y <- t(as.matrix(Y))
        rownames(Y) <- peaks[id1[1]]
        promoters_matrix_list[[i]] <- Y
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as.matrix(Y)
        target <- rownames(Y)
        regulators_matrix <- t(matrix)
        weight_list[[i]] <- sub.inferCSN(
          regulators_matrix = t(X),
          targets_matrix = t(matrix),
          target = target,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          k_folds = k_folds,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      }
    }
  } else {
    i <- NULL

    doParallel::registerDoParallel(cores = cores)
    if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weight_list <- foreach::foreach(
      i = seq_along(targets),
      .export = c("sub.inferCSN", "sparse.regression")
    ) %dopar% {
      p1 <- paste0(chr[i], ":", (starts[i] - 500), "-", starts[i])
      p2 <- paste0(chr[i], ":", (starts[i] - 250000), "-", (starts[i] + 250000))
      promoters <- cicero::find_overlapping_coordinates(peaks, p1)
      enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
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
        Y <- t(as.matrix(Y))
        rownames(Y) <- peaks[id1[1]]
        promoters_matrix_list[[i]] <- Y
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as.matrix(Y)
        target <- rownames(Y)
        regulators_matrix <- t(matrix)
        targets_matrix <- t(matrix)
        sub.inferCSN(
          regulators_matrix = X,
          targets_matrix = Y,
          target = target,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          k_folds = k_folds,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      }
    }
    doParallel::stopImplicitCluster()
  }

  weight_table <- purrr::list_rbind(weight_list)
  weight_table <- weight_table[order(abs(as.numeric(weight_table$weight)), decreasing = TRUE), ]
  # colnames(weight_table) <- c("Peak1", "Peak2", "weight")

  # variable selection
  if (is.null(high_corr_cutoff)) {
    high_corr_cutoff <- max(stats::quantile(weight_table$weight, 0.50), 0.001)
  }
  if (is.null(low_corr_cutoff)) {
    low_corr_cutoff <- min(0.001, stats::quantile(weight_table$weight, 0.25))
  }

  for (i in 1:length(weight_list)) {
    if (!is.null(weight_list[[i]])) {
      weight_table_h <- weight_list[[i]]
      coefficients <- weight_table_h[, 3]
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
          hic_index <- which(rownames(X) %in% weight_table_h[, 2][index1])
          other_index <- which(rownames(X) %in% weight_table_h[, 2][-index1])
          CPi_sub <- CPi[hic_index, other_index, drop = FALSE]
          flag_matrix <- matrix(0, nrow = nrow(CPi_sub), ncol = ncol(CPi_sub))
          flag_matrix[which(CPi_sub > 0.25)] <- 1
          correlated_index <- which(colSums(flag_matrix) > 0)
          if (!is.null(correlated_index)) {
            function_type[weight_table_h$Peak2 %in% rownames(X)[other_index[correlated_index]]] <- "high_corr"
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
          weight_table_h,
          function_type = function_type
        )
      )
    }
  }
  weight_table <- purrr::list_rbind(weight_list)
  weight_table$Starts <- as.numeric(weight_table$Starts)
  weight_table$Ends <- as.numeric(weight_table$Ends)
  weight_table$weight <- abs(as.numeric(weight_table$weight))

  return(weight_table)
}
