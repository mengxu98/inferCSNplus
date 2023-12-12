#' @param object An expression matrix, cells by genes
#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0" and "L0L2".
#' For high-dimensional and sparse data, such as single-cell sequencing data, "L0L2" is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param crossValidation Check whether cross validation is used.
#' @param nFolds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param kFolds The number of folds for sample split.
#' @param rThreshold rThreshold.
#' @param regulators Regulator genes.
#' @param targets Target genes.
#' @param maxSuppSize The number of non-zore coef, this value will affect the final performance.
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
    crossValidation = FALSE,
    seed = 1,
    nFolds = 10,
    kFolds = NULL,
    rThreshold = 0,
    regulators = NULL,
    targets = NULL,
    maxSuppSize = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
  if (verbose) message(paste("Running start for <", class(object)[1], ">."))
  matrix <- object
  # Check input parameters
  check.parameters(
    matrix = matrix,
    penalty = penalty,
    algorithm = algorithm,
    crossValidation = crossValidation,
    seed = seed,
    nFolds = nFolds,
    kFolds = kFolds,
    rThreshold = rThreshold,
    regulators = regulators,
    targets = targets,
    maxSuppSize = maxSuppSize,
    verbose = verbose,
    cores = cores
  )

  if (!is.null(regulators)) {
    regulators_matrix <- matrix[, intersect(colnames(matrix), regulators)]
  } else {
    regulators_matrix <- matrix
  }

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targetsMatrix <- matrix
  }
  targets <- colnames(targetsMatrix)
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

    weightDT <- purrr::map_dfr(targets, function(x) {
      if (verbose) pb$tick()
      sub.inferCSN(
        regulatorsMatrix = regulators_matrix,
        targetsMatrix = targetsMatrix,
        target = x,
        crossValidation = crossValidation,
        seed = seed,
        penalty = penalty,
        algorithm = algorithm,
        nFolds = nFolds,
        kFolds = kFolds,
        rThreshold = rThreshold,
        maxSuppSize = maxSuppSize,
        verbose = verbose
      )
    })
  } else {
    doParallel::registerDoParallel(cores = cores)
    if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")
    target <- NULL
    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(
      target = targets,
      .export = c("sub.inferCSN", "sparse.regression")
    ) %dopar% {
      sub.inferCSN(
        regulatorsMatrix = regulators_matrix,
        targetsMatrix = targetsMatrix,
        target = target,
        crossValidation = crossValidation,
        seed = seed,
        penalty = penalty,
        algorithm = algorithm,
        nFolds = nFolds,
        kFolds = kFolds,
        rThreshold = rThreshold,
        maxSuppSize = maxSuppSize,
        verbose = verbose
      )
    }
    weightDT <- data.table::rbindlist(weightDT)
    attr(weightDT, ".internal.selfref") <- NULL
    doParallel::stopImplicitCluster()
  }

  weightDT <- weightDT[order(abs(as.numeric(weightDT$weight)), decreasing = TRUE), ]
  if (verbose) message("Run done.")
  return(weightDT)
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
    crossValidation = FALSE,
    seed = 1,
    nFolds = 10,
    kFolds = NULL,
    rThreshold = 0,
    regulators = NULL,
    targets = NULL,
    maxSuppSize = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
  if (verbose) message(paste("Running start for <", class(object)[1], ">."))
  if (verbose) warning("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(object)
  inferCSN(
    matrix,
    penalty = penalty,
    algorithm = algorithm,
    crossValidation = crossValidation,
    seed = seed,
    nFolds = nFolds,
    kFolds = kFolds,
    rThreshold = rThreshold,
    regulators = regulators,
    targets = targets,
    maxSuppSize = maxSuppSize,
    verbose = verbose,
    cores = cores,
    ...
  )
}

#' @export
#' @method inferCSN Seurat
#'
#' @rdname inferCSN
#'
inferCSN.Seurat <- function(
    object,
    penalty = "L0",
    algorithm = "CD",
    crossValidation = FALSE,
    seed = 1,
    nFolds = 10,
    kFolds = NULL,
    rThreshold = 0,
    regulators = NULL,
    targets = NULL,
    maxSuppSize = NULL,
    verbose = FALSE,
    cores = 1,
    aggregate = TRUE,
    peakcalling = FALSE,
    macs2.path = NULL,
    fragments = NULL,
    k_neigh = 50,
    atacbinary = TRUE,
    max_overlap = 0.8,
    reduction.name = NULL,
    size_factor_normalize = FALSE,
    genome.info,
    early_stop = FALSE, # crossvaliation + sample?
    HC_cutoff = NULL,
    LC_cutoff = NULL,
    rescued = FALSE,
    ...) {
  if (verbose) message(paste("Running start for <", class(object)[1], ">."))

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

  # TODO: here indicate following analysis use the aggregated data,
  #   we should give a choice for users that they could decide whether aggregate data.
  if (aggregate) {
    # step 1. Aggregation
    if (verbose) {
      message("Generating aggregated data")
    }
    if ("aggregated_data" %in% names(Seurat::Misc(object))) {
      agg_data <- Seurat::Misc(object, slot = "aggregated_data")
    } else {
      agg_data <- aggregate.data(
        object,
        k_neigh = k_neigh,
        atacbinary = atacbinary,
        max_overlap = max_overlap,
        reduction.name = NULL,
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
      data_rna <- matrix(0, nrow = nrow(object@assays$RNA$counts), ncol = 1)
    }
    if ("ATAC" %in% names(object@assays)) {
      data_atac <- matrix(0, nrow = nrow(object@assays$ATAC@counts), ncol = 1)
    }
  }

  data_atac <- peaks.filter(data_atac)

  if (is.null(targets)) {
    # TODO: here code only could used to select targets for RNA, how to use for ATAC
    targets <- rownames(data_rna)
    targets <- na.omit(object@assays$RNA@meta.data$var.features)
  }

  ###
  if ("ATAC" %in% names(object@assays)) {
    rownames(data_atac) <- gsub("-", "_", rownames(data_atac))

    # Obtain candidate regions of focus markers
    genes <- lapply(genome.info$genes, function(x) strsplit(x, "[|]")[[1]][1])
    genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
    genes <- unlist(genes)
    genome.info$genes <- genes
    unik <- !duplicated(genes)
    genome.info <- genome.info[unik, ]

    if (is.null(targets)) {
      # TODO: here code only could used to select targets for RNA, how to use for ATAC
      targets <- genome.info$genes
    }

    targets <- lapply(targets, function(x) strsplit(x, "[.]")[[1]][1])
    targets <- unique(unlist(targets))
    targets <- genome.info$genes[which(genome.info$genes %in% targets)]

    genome.info.used <- genome.info[which(genome.info$genes %in% targets), ]
    chr <- genome.info.used$Chrom
    starts <- genome.info.used$Starts
    ends <- genome.info.used$Ends
  }

  if ("RNA" %in% names(agg_data)) {
    weight_network_rna <- inferCSN(
      matrix = t(data_rna),
      penalty = penalty,
      algorithm = algorithm,
      crossValidation = crossValidation,
      seed = seed,
      nFolds = nFolds,
      kFolds = kFolds,
      rThreshold = rThreshold,
      regulators = targets,
      targets = targets,
      maxSuppSize = maxSuppSize,
      verbose = verbose,
      cores = cores,
      ...
    )
  }

  if ("ATAC" %in% names(agg_data)) {
    weight_network_atac <- .inferCSN.atac(
      matrix = data_atac,
      penalty = penalty,
      algorithm = algorithm,
      crossValidation = crossValidation,
      seed = seed,
      nFolds = nFolds,
      kFolds = kFolds,
      rThreshold = rThreshold,
      regulators = targets,
      targets = targets,
      maxSuppSize = maxSuppSize,
      verbose = verbose,
      cores = cores,
      chr,
      starts,
      ...
    )
  }
  weightDT <- weight_network_atac
  weightDT <- weightDT[order(abs(as.numeric(weightDT$weight)), decreasing = TRUE), ]
  if (verbose) message("Run done.")
  return(weightDT)
}

.inferCSN.atac <- function(
    matrix,
    penalty = "L0",
    algorithm = "CD",
    crossValidation = FALSE,
    seed = 1,
    nFolds = 10,
    kFolds = NULL,
    rThreshold = 0,
    regulators = NULL,
    targets = NULL,
    maxSuppSize = NULL,
    verbose = FALSE,
    cores = 1,
    chr,
    starts,
    ...) {
  matrix <- as.matrix(matrix)

  if (!is.null(regulators)) {
    regulators_matrix <- matrix[, intersect(colnames(matrix), regulators)]
  } else {
    regulators_matrix <- matrix
  }

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targetsMatrix <- matrix
  }

  targets <- rownames(matrix)

  cores <- min(
    (parallel::detectCores(logical = FALSE) - 1), cores, length(targets)
  )
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
    TXs <- list()
    TYs <- list()
    weightDT <- purrr::map_dfr(
      seq_along(targets), function(i) {
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
          TXs[[i]] <- X
          Y <- matrix[id1, ]
          if (length(id1) > 1) {
            Y <- colSums(Y)
          }
          Y <- t(as.matrix(Y))
          rownames(Y) <- peaks[id1[1]] # Only use the first?
          TYs[[i]] <- Y
          if (ncol(X) == 1) {
            X <- t(X)
          }
          Y <- as.matrix(Y)
          target <- rownames(Y)
          regulators_matrix <- t(matrix)
          targetsMatrix <- t(matrix)
          sub.inferCSN(
            regulatorsMatrix = regulators_matrix,
            targetsMatrix = targetsMatrix,
            target = target,
            crossValidation = crossValidation,
            seed = seed,
            penalty = penalty,
            algorithm = algorithm,
            nFolds = nFolds,
            kFolds = kFolds,
            rThreshold = rThreshold,
            maxSuppSize = maxSuppSize,
            verbose = verbose
          )
        }
      }
    )
  } else {
    doParallel::registerDoParallel(cores = cores)
    if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    i <- NULL
    weightDT <- foreach::foreach(
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
        TXs[[i]] <- X
        Y <- matrix[id1, ]
        if (length(id1) > 1) {
          Y <- colSums(Y)
        }
        Y <- t(as.matrix(Y))
        rownames(Y) <- peaks[id1[1]] # Only use the first?
        TYs[[i]] <- Y
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as.matrix(Y)
        target <- rownames(Y)
        regulators_matrix <- t(matrix)
        targetsMatrix <- t(matrix)
        sub.inferCSN(
          regulatorsMatrix = regulators_matrix,
          targetsMatrix = targetsMatrix,
          target = target,
          crossValidation = crossValidation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          nFolds = nFolds,
          kFolds = kFolds,
          rThreshold = rThreshold,
          maxSuppSize = maxSuppSize,
          verbose = verbose
        )
      }
    }
    weightDT <- data.table::rbindlist(weightDT)
    attr(weightDT, ".internal.selfref") <- NULL
    doParallel::stopImplicitCluster()
  }

  weightDT <- weightDT[order(abs(as.numeric(weightDT$weight)), decreasing = TRUE), ]
  if (verbose) message("Run done.")
  return(weightDT)
}
