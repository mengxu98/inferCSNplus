.inferCSN.rna <- function(
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
    ...
) {
  warning("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(object)
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
  regulators <- colnames(regulators_matrix)

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
                                     total = length(regulators),
                                     clear = TRUE,
                                     width = 80)

    weightDT <- purrr::map_dfr(regulators, function(x) {
      if (verbose) pb$tick()
      sub.inferCSN(regulators_matrix = regulators_matrix,
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
                   verbose = verbose)
    })

  } else {
    doParallel::registerDoParallel(cores = cores)
    if(verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(target = targets,
                                 .export = c("sub.inferCSN", "sparse.regression")) %dopar% {
                                   sub.inferCSN(regulators_matrix = regulators_matrix,
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
    chr,
    starts,
    ...
) {
  warning("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(object)
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
  regulators <- colnames(regulators_matrix)

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
                                     total = length(regulators),
                                     clear = TRUE,
                                     width = 80)

    inferCSN_results <- list()
    TXs <- list()
    TYs <- list()

    weightDT <- purrr::map_dfr(seq_along(regulators), function(i) {
    # weightDT <- purrr::map_dfr(regulators, function(x) {
      if (verbose) pb$tick()
      # sub.inferCSN(regulators_matrix = regulators_matrix,
      #              targetsMatrix = targetsMatrix,
      #              target = x,
      #              crossValidation = crossValidation,
      #              seed = seed,
      #              penalty = penalty,
      #              algorithm = algorithm,
      #              nFolds = nFolds,
      #              kFolds = kFolds,
      #              rThreshold = rThreshold,
      #              maxSuppSize = maxSuppSize,
      #              verbose = verbose)
      if (verbose) {
        message("Inferring links for: ", focus_markers[i])
      }

      p1 <- paste(chr[i], ":", starts[i] - 500, "-", starts[i], sep = "")
      p2 <- paste(chr[i], ":", starts[i] - 250000, "-", starts[i] + 250000, sep = "")
      promoters <- cicero::find_overlapping_coordinates(peaks, p1)
      enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
      enhancers <- setdiff(enhancers, promoters)

      # if (verbose) {
      #   message("Inferring links for: ", focus_markers[i])
      # }

      ########## build data matrix of each gene used in gbm model ##########
      # if ("RNA" %in% names(agg_data)) {
      #   idx <- which(rownames(data_rna) == focus_markers[i])
      # } else {
      #   idx <- 1
      # }

      if (all((length(promoters) > 0 && length(enhancers) > 1))) {
        # if ((length(promoters) > 0 && length(enhancers) > 1) && length(idx) != 0) {
        id1 <- match(promoters, peaks)
        id1 <- id1[!is.na(id1)]
        id2 <- match(enhancers, peaks)
        id2 <- id2[!is.na(id2)]
        id2_new <- setdiff(id2, id1)
        X <- data_atac[id2_new, ]
        TXs[[i]] <- X
        Y <- data_atac[id1, ]
        if (length(id1) > 1) {
          Y <- colSums(Y)
        }
        Y <- t(as.matrix(Y))
        rownames(Y) <- peaks[id1[1]]
        TYs[[i]] <- Y
        # if ("RNA" %in% names(agg_data)) {
        #   Z <- data_rna[idx, ]
        #   Z <- t(as.matrix(Z))
        #   rownames(Z) <- focus_markers[i]
        # } else {
        #   Z <- Y
        # }
        Z <- Y
        flag <- 1
      } else {
        flag <- 0
        message("There are less than two peaks detected within 500 kb for ", focus_markers[i])
      }
      if (flag == 1) {
        X <- as.matrix(X)
        if (ncol(X) == 1) {
          X <- t(X)
        }
        Y <- as.matrix(Y)
        Z <- as.matrix(Z)
        rownames(Z) <- rownames(Y)
      }
      conns_h <- sub.inferCSN(regulators_matrix = regulators_matrix,
                   targetsMatrix = targetsMatrix,
                   target = focus_markers[i],
                   crossValidation = crossValidation,
                   seed = seed,
                   penalty = penalty,
                   algorithm = algorithm,
                   nFolds = nFolds,
                   kFolds = kFolds,
                   rThreshold = rThreshold,
                   maxSuppSize = maxSuppSize,
                   verbose = verbose)
      if (length(ncol(conns_h))) {
        inferCSN_results[[i]] <- conns_h
      }
    })

  } else {
    doParallel::registerDoParallel(cores = cores)
    if(verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(target = targets,
                                 .export = c("sub.inferCSN", "sparse.regression")) %dopar% {
                                   sub.inferCSN(regulators_matrix = regulators_matrix,
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
