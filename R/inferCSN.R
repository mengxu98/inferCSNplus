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
inferCSN.default <- function(object,
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
  if(verbose) message("Running start.")
  matrix <- object
  # Check input parameters
  check.parameters(matrix = matrix,
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
                   cores = cores)

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
#' @method inferCSN data.frame
#'
inferCSN.data.frame <- function(object,
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
  warning("Converting the class type of input data from <data.frame> to <matrix>.")
  matrix <- as.matrix(object)
  if (verbose) message("Running start.")

  # Check input parameters
  check.parameters(matrix = matrix,
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
                   cores = cores)

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
#' @method inferCSN Seurat
#'
inferCSN.Seurat <- function(object,
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
                            focus_markers = NULL,
                            early_stop = FALSE, # crossvaliation + sample?
                            HC_cutoff = NULL,
                            LC_cutoff = NULL,
                            rescued = FALSE,
                            seed = 123,
                             ...) {
  if (verbose) message("Running start.")

    # step 0. Peak calling

  if (("ATAC" %in% names(agg_data)) && peakcalling) {
    if (verbose) {
      message("Calling Peak")
    }
    object$cluster <- Idents(object)
    if (is.null(macs2.path)) {
      message("Please give the path to macs2!")
    }
    peaks <- CallPeaks(
      object = object,
      group.by = "cluster",
      macs2.path = macs2.path
    )
    if (is.null(fragments)) {
      message("Please input fragments!")
    }
    new_atac_data <- Signac:::FeatureMatrix(
      fragments = fragments, # fragments of original, we need to do aggregation
      features = peaks
    )
    object@assays$ATAC@counts <- new_atac_data
    if (verbose) {
      message("Peak calling finished")
    }
  }

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
      agg_data <- aggregate.data(object,
                                 k_neigh = k_neigh,
                                 atacbinary = atacbinary,
                                 max_overlap = max_overlap,
                                 reduction.name = NULL,
                                 size_factor_normalize = size_factor_normalize)
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

  if (is.null(focus_markers)) {
    focus_markers <- rownames(data_rna)
    focus_markers <- na.omit(object@assays$RNA@meta.data$var.features)
  }
  ###
  if ("ATAC" %in% names(object@assays)) {
    rownames(data_atac) <- gsub("-", "_", rownames(data_atac))
    peaks <- rownames(data_atac)

    # Obtain candidate regions of focus markers
    genes <- lapply(genome.info$genes, function(x) strsplit(x, "[|]")[[1]][1])
    genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
    genes <- unlist(genes)
    genome.info$genes <- genes
    unik <- !duplicated(genes)
    genome.info <- genome.info[unik, ]

    focus_markers <- lapply(focus_markers, function(x) strsplit(x, "[.]")[[1]][1])
    focus_markers <- unique(unlist(focus_markers))
    focus_markers <- genome.info$genes[which(genome.info$genes %in% focus_markers)]

    genome.info.used <- genome.info[which(genome.info$genes %in% focus_markers), ]
    chr <- genome.info.used$Chrom
    starts <- genome.info.used$Starts
    ends <- genome.info.used$Ends
  }

  if ("RNA" %in% names(agg_data)) {
    weight_network_rna <- .inferCSN.rna(object = t(data_rna),
                                        penalty = penalty,
                                        algorithm = algorithm,
                                        crossValidation = crossValidation,
                                        seed = 1,
                                        nFolds = nFolds,
                                        kFolds = kFolds,
                                        rThreshold = rThreshold,
                                        regulators = focus_markers,
                                        targets = targets,
                                        maxSuppSize = maxSuppSize,
                                        verbose = verbose,
                                        cores = cores)
  }

  if ("ATAC" %in% names(agg_data)) {
    weight_network_atac <- .inferCSN.atac(object = t(data_rna),
                                          penalty = penalty,
                                          algorithm = algorithm,
                                          crossValidation = crossValidation,
                                          seed = 1,
                                          nFolds = nFolds,
                                          kFolds = kFolds,
                                          rThreshold = rThreshold,
                                          regulators = focus_markers,
                                          targets = targets,
                                          maxSuppSize = maxSuppSize,
                                          verbose = verbose,
                                          cores = cores,
                                          chr,
                                          starts)
  }

  for (i in 1:length(focus_markers)) {



      ### whether use early stop rule
      if (infer_method == "SRM") {
        if (early_stop) {
          if (length(size_factor) / 5 < 100) {
            message("The number of cells is too small to split!")
          }
          set.seed(seed)
          cv_idx <- sample(1:5, size = length(size_factor), replace = T)
          test_idx <- which(cv_idx == 1)
          validation_idx <- which(cv_idx == 2)
          x_train <- as.matrix(X[, -c(test_idx, validation_idx)])
          x_test <- as.matrix(X[, test_idx])
          x_validation <- as.matrix(X[, validation_idx])
          y_train <- Y[-c(test_idx, validation_idx)]
          y_test <- Y[test_idx]
          y_validation <- Y[validation_idx]

          dtrain <- xgb.DMatrix(data = t(x_train), label = as.numeric(y_train))
          dtest <- xgb.DMatrix(data = t(x_test), label = as.numeric(y_test))
          dvalidation <- xgb.DMatrix(data = t(x_validation), label = as.numeric(y_validation))
          watchlist1 <- list(train = dtrain, test = dvalidation)
          watchlist2 <- list(train = dtrain, test = dtest)
          xgb_v <- xgb.train(
            params = params,
            data = dtrain,
            watchlist = watchlist1,
            nrounds = 100,
            nthread = nthread,
            objective = "reg:linear",
            verbose = 0
          )
          cv1 <- xgb_v$evaluation_log
          rmse_d <- cv1$test_rmse - cv1$train_rmse
          rmse_dd <- abs(rmse_d[2:100] - rmse_d[1:99]) / rmse_d[1]
          stop_index <- which(rmse_dd == min(rmse_dd))
          # train final model
          xgb.fit.final <- xgboost(
            params = params,
            data = t(X),
            label = as.numeric(Z),
            nrounds = stop_index[1],
            nthread = nthread,
            objective = "reg:squarederror",
            verbose = 0
          )
        } else {
          # train final model
          xgb.fit.final <- xgboost(
            params = params,
            data = t(X),
            label = as.numeric(Z),
            nrounds = 100,
            nthread = nthread,
            objective = "reg:squarederror",
            verbose = 0
          )
        }

        # create importance matrix
        tryCatch(
          {
            importance_matrix <- xgb.importance(model = xgb.fit.final)
            Imp_peak <- importance_matrix$Feature
            Imp_peak <- as.vector(Imp_peak)
            Imp_value <- importance_matrix$Gain
            Imp_peak_h <- Imp_peak
            Imp_value_h <- Imp_value
            conns_h <- list()
            conns_h$Peak1 <- as.character(rownames(Y))
            conns_h$Peak2 <- as.character(Imp_peak_h)
            conns_h$Importance <- Imp_value_h
            conns_h <- as.data.frame(conns_h)
            colnames(conns_h) <- c("Peak1", "Peak2", "Importance")
          },
          error = function(e) {
          }
        )
      }



    }
  }
  conns <- do.call(rbind, inferCSN_results)
  ######################### variable selection ###########
  if (is.null(HC_cutoff)) {
    HC_cutoff <- max(stats::quantile(conns$Importance, 0.50), 0.001)
  }
  if (is.null(LC_cutoff)) {
    LC_cutoff <- min(0.001, stats::quantile(conns$Importance, 0.25))
  }
  for (i in 1:length(inferCSN_results)) {
    if (!is.null(inferCSN_results[[i]])) {
      conns_h <- inferCSN_results[[i]]
      Imp_value <- conns_h$Importance
      index1 <- which(Imp_value > HC_cutoff) # HC
      index2 <- intersect(which(Imp_value > LC_cutoff), which(Imp_value <= HC_cutoff)) # MC
      index3 <- which(Imp_value <= LC_cutoff) # LC

      #### do data frame:gene, starts, end, peak1,peak2,importance,function_type
      function_type <- rep(NA, length(Imp_value))
      function_type[index1] <- "HC"
      function_type[index2] <- "MC"
      function_type[index3] <- "LC"
      # rescue highly correlated CREs
      if (rescued) {
        if (i <= length(TXs)) {
          X <- TXs[[i]]
          Y <- TYs[[i]]
          CPi <- abs(cor(t(X)))
          for (p in 1:nrow(CPi)) {
            CPi[p, p] <- 0
          }
          # focus on HC rows
          hic_index <- which(rownames(X) %in% conns_h$Peak2[index1])
          other_index <- which(rownames(X) %in% conns_h$Peak2[-index1])
          CPi_sub <- CPi[hic_index, other_index, drop = FALSE]
          flag_matrix <- matrix(0, nrow = nrow(CPi_sub), ncol = ncol(CPi_sub))
          flag_matrix[which(CPi_sub > 0.25)] <- 1
          correlated_index <- which(colSums(flag_matrix) > 0)
          if (!is.null(correlated_index)) {
            function_type[conns_h$Peak2 %in% rownames(X)[other_index[correlated_index]]] <- "HC"
          }
        }
      }
      inferCSN_results[[i]] <- cbind(data.frame(gene = focus_markers[i],
                                                Chr = Chr[i],
                                                Starts = Starts[i],
                                                Ends = Ends[i]),
                                     cbind(conns_h,
                                           function_type = function_type))
    }
  }
  inferCSN_Result_all <- do.call(rbind, inferCSN_results)
  inferCSN_Result_all$Starts <- as.numeric(inferCSN_Result_all$Starts)
  inferCSN_Result_all$Ends <- as.numeric(inferCSN_Result_all$Ends)
  inferCSN_Result_all$Importance <- as.numeric(inferCSN_Result_all$Importance)
  # save result
  Seurat::Misc(object, slot = "direct.net") <- inferCSN_Result_all
  # return(object)

  ###

  # Check input parameters
  check.parameters(matrix = matrix,
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
                   cores = cores)

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
