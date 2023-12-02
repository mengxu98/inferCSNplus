#' run inferCSN on Seurat object
#'
#' @param object Seurat object.
#' @param peakcalling  call peak
#' @param macs2.path  path to macs2
#' @param fragments  fragments file
#' @param k_neigh Number of cells to aggregate per group.
#' @param atacbinary Logical, should accessibility values be binarized
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param reduction.name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, whether need to do size normalization
#' @param genome.info the TSS information of genome, e.g. hg19, hg38
#' @param focus_markers the focused genes
#' @param params the list of parameters used in Xgboost
#' @param nthread  the number of threads can be manually specified in Xgboost trainning stage, default is 2
#' @param early_stop Logical, whether use early stop rule on validation data to reduce overfitting
#' @param HC_cutoff the threshold of high functional CREs
#' @param LC_cutoff the threshold of low functional CREs
#' @param rescued Logical, whether to rescue highly correlated CREs
#' @param seed Random seed
#' @param verbose Logical, should warning and info messages be printed
#' @import Signac
#' @import Seurat
#' @import Matrix
#' @import xgboost
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom cicero find_overlapping_coordinates
#' @return a Seurat object with new links assay.
#' @export
inferCSNplus <- function(object,
                         peakcalling = FALSE,
                         macs2.path = NULL,
                         fragments = NULL,
                         k_neigh = 50,
                         atacbinary = TRUE,
                         max_overlap = 0.8,
                         reduction.name = NULL,
                         size_factor_normalize = FALSE,
                         genome.info,
                         focus_markers,
                         params = NULL,
                         nthread = 2,
                         early_stop = FALSE,
                         HC_cutoff = NULL,
                         LC_cutoff = NULL,
                         rescued = FALSE,
                         seed = 123,
                         verbose = TRUE,
                         infer_method = "SRM") {
  # step 0. Peak calling

  if (peakcalling) {
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
    new_atac_data <- FeatureMatrix(
      fragments = fragments, # fragments of original, we need to do aggregation
      features = peaks
    )
    object@assays$ATAC@counts <- new_atac_data
    if (verbose) {
      message("Peak calling finished")
    }
  }

  # step 1. Aggregation
  if (verbose) {
    message("Generating aggregated data")
  }
  if ("aggregated.data" %in% names(Seurat::Misc(object))) {
    agg_data <- Seurat::Misc(object, slot = "aggregated.data")
  } else {
    agg_data <- aggregate.data(object,
                               k_neigh = k_neigh,
                               atacbinary = atacbinary,
                               max_overlap = max_overlap,
                               reduction.name = NULL,
                               size_factor_normalize = size_factor_normalize)
    Seurat::Misc(object, slot = "aggregated.data") <- agg_data
  }

  # step2.  inferCSNplus
  # library(cicero)
  # library(Matrix)
  # library(data.table)
  ## import gbm package
  # library(xgboost)      # a faster implementation of gbm
  # library(ggplot2)      # model visualization
  options(stringsAsFactors = FALSE)

  ########## Implement gbm for each gene  ##########
  # parameter list for xgboost
  if (is.null(params)) {
    params <- list(
      eta = 0.3,
      max_depth = 6,
      min_child_weight = 1,
      subsample = 1,
      colsample_bytree = 1,
      lambda = 1
    )
  }

  # TODO: here indicate following analysis use the aggregated data,
  #   we should give a choice for users that they could decide whether aggregate data.
  if ("RNA" %in% names(agg_data)) {
    data_rna <- as.matrix(agg_data$RNA)

    genes <- rownames(data_rna)
    # genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
    # genes <- unlist(genes)
    # genes <- toupper(genes)
    rownames(data_rna) <- genes
    data_rna <- data_rna[unique(genes), ]
  }

  if ("ATAC" %in% names(agg_data)) {
    data_atac <- as.matrix(agg_data$ATAC)
    rownames(data_atac) <- gsub("-", "_", rownames(data_atac))
    peaks <- rownames(data_atac)

    ########## Obtain candidate regions of focus markers ##########
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
    Chr <- genome.info.used$Chrom
    Starts <- genome.info.used$Starts
    Ends <- genome.info.used$Ends
  }

  inferCSN_results <- list()
  TXs <- list()
  TYs <- list()
  for (i in 1:length(focus_markers)) {
    if (verbose) {
      message("Inferring links for: ", focus_markers[i])
    }

    p1 <- paste(Chr[i], ":", Starts[i] - 500, "-", Starts[i], sep = "")
    p2 <- paste(Chr[i], ":", Starts[i] - 250000, "-", Starts[i] + 250000, sep = "")
    promoters <- cicero::find_overlapping_coordinates(peaks, p1)
    enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
    enhancers <- setdiff(enhancers, promoters)

    ########## build data matrix of each gene used in gbm model ##########
    if ("RNA" %in% names(agg_data)) {
      idx <- which(rownames(data_rna) == focus_markers[i])
    } else {
      idx <- 1
    }

    if ((length(promoters) > 0 && length(enhancers) > 1) && length(idx) != 0) {
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
      if ("RNA" %in% names(agg_data)) {
        Z <- data_rna[idx, ]
        Z <- t(as.matrix(Z))
        rownames(Z) <- focus_markers[i]
      } else {
        Z <- Y
      }
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
      ### whether use early stop rule
      if (infer_method == "SRM") {

        print("Using SRM.")
      } else {
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


      if (length(ncol(conns_h))) {
        inferCSN_results[[i]] <- conns_h
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
  return(object)
}
