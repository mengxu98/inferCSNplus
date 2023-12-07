.inferCSN.rna <- function(
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
    ...
) {
  warning("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(matrix)
  if (verbose) message("Running start.")

  # Check input parameters
  # check.parameters(matrix = matrix,
  #                  penalty = penalty,
  #                  algorithm = algorithm,
  #                  crossValidation = crossValidation,
  #                  seed = seed,
  #                  nFolds = nFolds,
  #                  kFolds = kFolds,
  #                  rThreshold = rThreshold,
  #                  regulators = regulators,
  #                  targets = targets,
  #                  maxSuppSize = maxSuppSize,
  #                  verbose = verbose,
  #                  cores = cores)

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
    if(verbose) message("Using 1 core.")
    # Format progress information
    format <- "Running [:bar] :percent, No.:current of :total genes, :elapsed."
    pb <- progress::progress_bar$new(format = format,
                                     total = length(targets),
                                     clear = TRUE,
                                     width = 80)

    weightDT <- c()
    for (i in seq_along(targets)) {
      x <- targets[i]
      if (verbose) pb$tick()
      weightDT <- rbind(weightDT, sub.inferCSN(regulatorsMatrix = regulators_matrix,
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
                                               verbose = verbose))

    }

    # weightDT <- purrr::map_dfr(targets, function(x) {
    #   if (verbose) pb$tick()
    #   sub.inferCSN(regulatorsMatrix = regulators_matrix,
    #                targetsMatrix = targetsMatrix,
    #                target = x,
    #                crossValidation = crossValidation,
    #                seed = seed,
    #                penalty = penalty,
    #                algorithm = algorithm,
    #                nFolds = nFolds,
    #                kFolds = kFolds,
    #                rThreshold = rThreshold,
    #                maxSuppSize = maxSuppSize,
    #                verbose = verbose)
    # })

  } else {
    doParallel::registerDoParallel(cores = cores)
    if(verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(target = targets,
                                 .export = c("sub.inferCSN", "sparse.regression")) %dopar% {
                                   sub.inferCSN(regulatorsMatrix = regulators_matrix,
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
                                                verbose = verbose)
                                 }
    weightDT <- data.table::rbindlist(weightDT)
    attr(weightDT, ".internal.selfref") <- NULL
    doParallel::stopImplicitCluster()
  }

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
    ...
) {
  matrix <- as.matrix(matrix)
  if (verbose) message("Running start.")

  # Check input parameters
  # check.parameters(matrix = matrix,
  #                  penalty = penalty,
  #                  algorithm = algorithm,
  #                  crossValidation = crossValidation,
  #                  seed = seed,
  #                  nFolds = nFolds,
  #                  kFolds = kFolds,
  #                  rThreshold = rThreshold,
  #                  regulators = regulators,
  #                  targets = targets,
  #                  maxSuppSize = maxSuppSize,
  #                  verbose = verbose,
  #                  cores = cores)

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
  # targets <- colnames(targetsMatrix)
  # rm(matrix)

  targets <- rownames(matrix)

  input_data_sum <- list()
  for (i in seq_along(targets)) {
    input_data <- list()
    p1 <- paste0(chr[i], ":", starts[i] - 500, "-", starts[i])
    p2 <- paste0(chr[i], ":", starts[i] - 250000, "-", starts[i] + 250000)
    promoters <- cicero::find_overlapping_coordinates(peaks, p1)
    enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
    enhancers <- setdiff(enhancers, promoters)

    if (length(promoters) > 0 && length(enhancers) > 1) {
      id1 <- na.omit(match(promoters, peaks))
      # id1 <- id1[!is.na(id1)]
      id2 <- na.omit(match(enhancers, peaks))
      # id2 <- id2[!is.na(id2)]
      X <- matrix[setdiff(id2, id1), ]

      TXs[[i]] <- X
      Y <- matrix[id1, ]
      if (length(id1) > 1) {
        Y <- colSums(Y)
      }
      Y <- t(as.matrix(Y))
      rownames(Y) <- peaks[id1[1]] # Only use the first?
      TYs[[i]] <- Y
      # X <- as.matrix(X)
      if (ncol(X) == 1) {
        X <- t(X)
      }
      Y <- as.matrix(Y)
      target <- rownames(Y)
      regulators_matrix <- t(matrix)
      targetsMatrix <- t(matrix)
    } else {
      next
    }
  }

  cores <- min((parallel::detectCores(logical = FALSE) - 1), cores, length(targets))
  if (cores == 1) {
    if(verbose) message("Using 1 core.")
    # Format progress information
    format <- "Running [:bar] :percent, No.:current of :total peaks, :elapsed."
    pb <- progress::progress_bar$new(format = format,
                                     total = length(targets),
                                     clear = TRUE,
                                     width = 80)

    peaks <- rownames(matrix)
    TXs <- list()
    TYs <- list()
    weightDT <- c()
    # for (i in seq_along(targets)) {
    #   if (verbose) pb$tick()
    #
    #   p1 <- paste0(chr[i], ":", starts[i] - 500, "-", starts[i])
    #   p2 <- paste0(chr[i], ":", starts[i] - 250000, "-", starts[i] + 250000)
    #   promoters <- cicero::find_overlapping_coordinates(peaks, p1)
    #   enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
    #   enhancers <- setdiff(enhancers, promoters)
    #
    #   if (length(promoters) > 0 && length(enhancers) > 1) {
    #     id1 <- na.omit(match(promoters, peaks))
    #     # id1 <- id1[!is.na(id1)]
    #     id2 <- na.omit(match(enhancers, peaks))
    #     # id2 <- id2[!is.na(id2)]
    #     X <- matrix[setdiff(id2, id1), ]
    #     TXs[[i]] <- X
    #     Y <- matrix[id1, ]
    #     if (length(id1) > 1) {
    #       Y <- colSums(Y)
    #     }
    #     Y <- t(as.matrix(Y))
    #     rownames(Y) <- peaks[id1[1]] # Only use the first?
    #     TYs[[i]] <- Y
    #     # X <- as.matrix(X)
    #     if (ncol(X) == 1) {
    #       X <- t(X)
    #     }
    #     Y <- as.matrix(Y)
    #     target <- rownames(Y)
    #     regulators_matrix <- t(matrix)
    #     targetsMatrix <- t(matrix)
    #     weightDT <- rbind(
    #       weightDT,
    #       sub.inferCSN(
    #         regulatorsMatrix = regulators_matrix,
    #         targetsMatrix = targetsMatrix,
    #         target = target,
    #         crossValidation = crossValidation,
    #         seed = seed,
    #         penalty = penalty,
    #         algorithm = algorithm,
    #         nFolds = nFolds,
    #         kFolds = kFolds,
    #         rThreshold = rThreshold,
    #         maxSuppSize = maxSuppSize,
    #         verbose = verbose
    #       )
    #     )
    #   } else {
    #     next
    #   }
    # }


    weightDT <- purrr::map_dfr(
      seq_along(targets), function(i) {
        if (verbose) pb$tick()

        p1 <- paste0(chr[i], ":", starts[i] - 500, "-", starts[i])
        p2 <- paste0(chr[i], ":", starts[i] - 250000, "-", starts[i] + 250000)
        promoters <- cicero::find_overlapping_coordinates(peaks, p1)
        enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
        enhancers <- setdiff(enhancers, promoters)

        if (length(promoters) > 0 && length(enhancers) > 1) {
          id1 <- na.omit(match(promoters, peaks))
          # id1 <- id1[!is.na(id1)]
          id2 <- na.omit(match(enhancers, peaks))
          # id2 <- id2[!is.na(id2)]
          X <- matrix[setdiff(id2, id1), ]
          TXs[[i]] <- X
          Y <- matrix[id1, ]
          if (length(id1) > 1) {
            Y <- colSums(Y)
          }
          Y <- t(as.matrix(Y))
          rownames(Y) <- peaks[id1[1]] # Only use the first?
          TYs[[i]] <- Y
          # X <- as.matrix(X)
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
    if(verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(
      i = seq_along(targets),
      .export = c("sub.inferCSN", "sparse.regression")) %dopar% {
        p1 <- paste0(chr[i], ":", starts[i] - 500, "-", starts[i])
        p2 <- paste0(chr[i], ":", starts[i] - 250000, "-", starts[i] + 250000)
        promoters <- cicero::find_overlapping_coordinates(peaks, p1)
        enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
        enhancers <- setdiff(enhancers, promoters)

        if (length(promoters) > 0 && length(enhancers) > 1) {
          id1 <- na.omit(match(promoters, peaks))
          # id1 <- id1[!is.na(id1)]
          id2 <- na.omit(match(enhancers, peaks))
          # id2 <- id2[!is.na(id2)]
          X <- matrix[setdiff(id2, id1), ]
          TXs[[i]] <- X
          Y <- matrix[id1, ]
          if (length(id1) > 1) {
            Y <- colSums(Y)
          }
          Y <- t(as.matrix(Y))
          rownames(Y) <- peaks[id1[1]] # Only use the first?
          TYs[[i]] <- Y
          # X <- as.matrix(X)
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
