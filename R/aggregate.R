#' Create aggregated data for inferCSN
#'
#' @description
#' Function to generate aggregated inputs for inferCSN. \code{aggregating_data}
#'
#' @inheritParams inferCSN
#'
#' @return Seurat object with aggregated data
#' @export
aggregating_data <- function(
    object,
    k_neigh = 50,
    atacbinary = TRUE,
    max_overlap = 0.8,
    reduction_name = NULL,
    size_factor_normalize = TRUE,
    seed = 123,
    verbose = FALSE) {
  if (verbose == 2) {
    message("---Setting 'aggregate = TRUE', and there is no aggregated data in seurat object.")
    message("---Aggregating data, and using 'names(Seurat::Misc(object))' to check it.")
  }
  if (!is.null(reduction_name)) {
    cell_coord <- object@reductions[[reduction_name]]
  } else {
    if (all(c("RNA", "ATAC") %in% names(object@assays))) {
      cell_coord <- object@reductions$wnn.umap@cell.embeddings
    } else {
      if ("umap" %in% names(object@reductions)) {
        cell_coord <- object@reductions$umap@cell.embeddings
      }
    }
  }

  cluster <- as.character(Seurat::Idents(object))
  uniqe_cluster <- unique(cluster)
  if ("RNA" %in% names(object@assays)) {
    rna_new <- matrix(0, nrow = nrow(object@assays$RNA$counts), ncol = 1)
  }
  if ("ATAC" %in% names(object@assays)) {
    atac_new <- matrix(0, nrow = nrow(object@assays$ATAC@counts), ncol = 1)
  }

  cell_sample <- matrix(0, nrow = 1, ncol = k_neigh)

  for (i in seq_along(uniqe_cluster)) {
    if (verbose) {
      message("Aggregating data for cluster: ", uniqe_cluster[i])
    }
    subobject <- subset(object, idents = uniqe_cluster[i])
    sub_index <- which(cluster %in% uniqe_cluster[i])
    cell_coord_i <- as.data.frame(cell_coord[sub_index, ])
    sub_aggregated_data <- generate_aggregated_data(
      subobject,
      cell_coord_i,
      k_neigh,
      atacbinary,
      max_overlap,
      seed,
      verbose
    )

    sub_cell_sample <- sub_aggregated_data$cell_sample
    if ("RNA" %in% names(object@assays)) {
      rna_new <- cbind(rna_new, sub_aggregated_data$RNA)
    }
    if ("ATAC" %in% names(object@assays)) {
      atac_new <- cbind(atac_new, sub_aggregated_data$ATAC)
    }
    if (ncol(sub_cell_sample) < k_neigh) {
      sub_cell_sample_new <- as.matrix(sub_cell_sample)
      sub_cell_sample_new <- cbind(
        sub_cell_sample_new,
        matrix(0, nrow = 1, ncol = k_neigh - ncol(sub_cell_sample_new))
      )
    } else {
      sub_cell_sample_new <- apply(sub_cell_sample, 2, function(x) {
        sub_index[x] # for each column return original index
      })
      sub_cell_sample_new <- as.data.frame(sub_cell_sample_new)
      sub_cell_sample_new <- as.matrix(sub_cell_sample_new)
    }
    cell_sample <- rbind(cell_sample, sub_cell_sample_new)
  }
  if ("RNA" %in% names(object@assays)) {
    rna_new <- rna_new[, -1]
  }
  if ("ATAC" %in% names(object@assays)) {
    atac_new <- atac_new[, -1]
    cell_sample <- cell_sample[-1, ]
  }

  # Normalization
  if (size_factor_normalize) {
    if ("RNA" %in% names(object@assays)) {
      rna_new <- t(t(log(rna_new + 1)) / estimate.size.factors(rna_new))
    }
    if ("ATAC" %in% names(object@assays)) {
      atac_new <- t(t(log(atac_new + 1)) / estimate.size.factors(atac_new))
    }
  }
  new_data <- list()
  if ("RNA" %in% names(object@assays)) {
    new_data$RNA <- rna_new
  }
  if ("ATAC" %in% names(object@assays)) {
    new_data$ATAC <- atac_new
  }

  new_data$cell_sample <- cell_sample
  return(new_data)
}

#' Create aggregated data for a certain cluster
#'
#' Function to generate aggregated inputs of a cetrain cluster. \code{generate.aggregated.data}
#' takes as input sparse data. This function will aggregate binary accessibility scores (or gene expression)
#' per cell cluster, if they do not overlap any existing cluster with more than 50% cells.
#'
#' @param object Seurat object.
#' @param cell_coord similarity matrix or dimiension reductions.
#' @param k_neigh Number of cells to aggregate per cluster.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary
#' @param max_overlap The maximum overlapping ratio of two clusters.
#' @param seed Random seed
#' @param verbose Logical, should warning and info messages be printed?
#'
#' @return Aggregated data.
#' @export
generate_aggregated_data <- function(
    object,
    cell_coord,
    k_neigh = 50,
    atacbinary = TRUE,
    max_overlap = 0.8,
    seed = 1,
    verbose = TRUE) {
  if (nrow(cell_coord) > k_neigh) {
    # Create a k-nearest neighbors map
    nn_map <- as.data.frame(
      FNN::knn.index(cell_coord, k = (k_neigh - 1))
    )
    row.names(nn_map) <- row.names(cell_coord)
    nn_map$agg_cell <- 1:nrow(nn_map)
    good_choices <- 1:nrow(nn_map)

    # Sample cells randomly
    set.seed(seed)
    choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]

    it <- 0
    # Slow (contain calculating of overlapping between cell clusters)
    while (length(good_choices) > 0 & it < nrow(cell_coord) / ((1 - max_overlap) * k_neigh)) {
      it <- it + 1
      choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]

      # calculate overlapping between cell clusters
      combs <- data.frame(1:(nrow(cell_sample) - 1), nrow(cell_sample))
      shared <- apply(combs, 1, function(x) { # Slow
        (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x, ]))))
      })
      if (max(shared) < max_overlap * k_neigh) {
        chosen <- new_chosen
      }
    }

    # Aggregating both scRNA-seq and scATAC-seq counts of cells within one cluster
    if ("RNA" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$RNA$counts)

      rna_mask <- sapply(
        seq_len(nrow(cell_sample)), function(x) {
          seq_len(ncol(rna_old)) %in% cell_sample[x, , drop = FALSE]
        }
      )

      rna_mask <- Matrix::Matrix(rna_mask)
      rna_new <- rna_old %*% rna_mask
      rna_new <- as.matrix(rna_new)
    }

    if ("ATAC" %in% names(object@assays)) {
      atac_old <- object@assays$ATAC@counts

      if (atacbinary) {
        atac_old <- atac_old > 0
      }

      atac_mask <- sapply(
        seq_len(nrow(cell_sample)), function(x) {
          seq_len(ncol(atac_old)) %in% cell_sample[x, , drop = FALSE]
        }
      )

      atac_mask <- Matrix::Matrix(atac_mask)
      atac_new <- atac_old %*% atac_mask
      atac_new <- as.matrix(atac_new)
    }
  } else {
    if ("RNA" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$RNA$counts)
      rna_new <- rowSums(rna_old)
      rna_new <- as.matrix(rna_new)
    }
    if ("ATAC" %in% names(object@assays)) {
      atac_old <- object@assays$ATAC@counts

      if (atacbinary) atac_old <- atac_old > 0

      atac_new <- rowSums(atac_old)
      atac_new <- as.matrix(atac_new)
    }

    cell_sample <- as.data.frame(t(matrix(seq(from = 1, to = nrow(cell_coord)))))
  }
  new_data <- list()
  if ("RNA" %in% names(object@assays)) {
    new_data$RNA <- rna_new
  }
  if ("ATAC" %in% names(object@assays)) {
    new_data$ATAC <- atac_new
  }
  new_data$cell_sample <- cell_sample

  return(new_data)
}

#' Function to calculate the size factor for the single-cell data
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches.
#' It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#' @importFrom stats median
#' @import slam
#' @import Matrix
#' @export
estimate.size.factors <- function(
    counts,
    locfunc = median,
    round_exprs = TRUE,
    method = "mean-geometric-mean-total") {
  if (is.sparse(counts)[1]) {
    estimate.size.factors.sparse(
      counts,
      locfunc = locfunc,
      round_exprs = round_exprs,
      method = method
    )
  } else {
    estimate.size.factors.dense(
      counts,
      locfunc = locfunc,
      round_exprs = round_exprs,
      method = method
    )
  }
}

#' Convert a slam matrix to a sparseMatrix
#' @param slam_matrix A slam matrix
#' @export
as_sparse <- function(slam_matrix) {
  retVal <- Matrix::sparseMatrix(
    i = slam_matrix[["i"]],
    j = slam_matrix[["j"]],
    x = slam_matrix[["v"]],
    dims = c(
      slam_matrix[["nrow"]],
      slam_matrix[["ncol"]]
    )
  )
  if (!is.null(slam_matrix[["dimnames"]])) {
    dimnames(retVal) <- slam_matrix[["dimnames"]]
  }
  return(retVal)
}

#' Convert a sparseMatrix from Matrix package to a slam matrix
#' @param sp_mat The matrix for the aggregated single cell data
#' @import slam
#' @export
as.slam.matrix <- function(sp_mat) {
  sp <- Matrix::summary(sp_mat)
  slam::simple_triplet_matrix(
    sp[, "i"],
    sp[, "j"],
    sp[, "x"],
    ncol = ncol(sp_mat),
    nrow = nrow(sp_mat),
    dimnames = dimnames(sp_mat))
}

#' is.sparse
#'
#' @param x Data
#' @export
is.sparse <- function(x) {
  class(x) %in% c("dgcountsatrix", "dgTMatrix")
}

#' Estimate size factors for each column, given a sparseMatrix from the Matrix package
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts.
#' @param locfunc The location function used to find the representive value.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded.
#' @param method A character to specify the size factor calculation appraoches.
#' It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
#' @importFrom stats median
#'
#' @export
estimate.size.factors.sparse <- function(
    counts,
    locfunc = median,
    round_exprs = TRUE,
    method = "mean-geometric-mean-total") {
  if (round_exprs) {
    counts <- round(counts)
  }
  counts <- as.slam.matrix(counts)

  if (method == "weighted-median") {
    log_medians <- slam::rowapply_simple_triplet_matrix(
      counts, function(cell_expr) {
        log(locfunc(cell_expr))
      }
    )

    weights <- slam::rowapply_simple_triplet_matrix(
      counts, function(cell_expr) {
        num_pos <- sum(cell_expr > 0)
        num_pos / length(cell_expr)
      }
    )

    sfs <- slam::colapply_simple_triplet_matrix(
      counts, function(cnts) {
        norm_cnts <- weights * (log(cnts) - log_medians)
        norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
        norm_cnts <- norm_cnts[is.finite(norm_cnts)]
        exp(mean(norm_cnts))
      }
    )
  } else if (method == "median-geometric-mean") {
    log_geo_means <- slam::rowapply_simple_triplet_matrix(
      counts, function(x) {
        mean(log(counts))
      }
    )

    sfs <- slam::colapply_simple_triplet_matrix(
      counts, function(cnts) {
        norm_cnts <- log(cnts) - log_geo_means
        norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
        norm_cnts <- norm_cnts[is.finite(norm_cnts)]
        exp(locfunc(norm_cnts))
      }
    )
  } else if (method == "median") {
    stop("Error: method 'median' not yet supported for sparse matrices")
  } else if (method == "mode") {
    stop("Error: method 'mode' not yet supported for sparse matrices")
  } else if (method == "geometric-mean-total") {
    cell_total <- slam::col_sums(counts)
    sfs <- log(cell_total) / mean(log(cell_total))
  } else if (method == "mean-geometric-mean-total") {
    cell_total <- slam::col_sums(counts)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1

  return(sfs)
}

#' Estimate size factors dense matrix
#'
#' @inheritParams estimate.size.factors.sparse
#'
#' @importFrom stats median
#'
#' @export
estimate.size.factors.dense <- function(
    counts,
    locfunc = median,
    round_exprs = TRUE,
    method = "mean-geometric-mean-total") {
  counts <- counts
  if (round_exprs) {
    counts <- round(counts)
  }
  if (method == "weighted-median") {
    log_medians <- apply(
      counts, 1, function(cell_expr) {
        log(locfunc(cell_expr))
      }
    )

    weights <- apply(
      counts, 1, function(cell_expr) {
        num_pos <- sum(cell_expr > 0)
        num_pos / length(cell_expr)
      }
    )

    sfs <- apply(
      counts, 2, function(cnts) {
        norm_cnts <- weights * (log(cnts) - log_medians)
        norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
        norm_cnts <- norm_cnts[is.finite(norm_cnts)]
        exp(mean(norm_cnts))
      }
    )
  } else if (method == "median-geometric-mean") {
    log_geo_means <- rowMeans(log(counts))

    sfs <- apply(
      counts, 2, function(cnts) {
        norm_cnts <- log(cnts) - log_geo_means
        norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
        norm_cnts <- norm_cnts[is.finite(norm_cnts)]
        exp(locfunc(norm_cnts))
      }
    )
  } else if (method == "median") {
    row_median <- apply(counts, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(counts) - row_median), 2, median)
  } else if (method == "mode") {
    sfs <- estimate.t(counts)
  } else if (method == "geometric-mean-total") {
    cell_total <- apply(counts, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  } else if (method == "mean-geometric-mean-total") {
    cell_total <- apply(counts, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  return(sfs)
}

#' Find the most commonly occuring relative expression value in each cell
#'
#' @description Converting relative expression values to mRNA copies per cell requires
#' knowing the most commonly occuring relative expression value in each cell
#' This value typically corresponds to an RPC value of 1. This function
#' finds the most commonly occuring (log-transformed) relative expression value
#' for each column in the provided expression matrix.
#'
#' @details This function estimates the most abundant relative expression value
#' (t^*) using a gaussian kernel density function. It can also optionally
#' output the t^* based on a two gaussian mixture model
#' based on the smsn.mixture from mixsmsn package
#'
#' @param relative_expr_matrix a matrix of relative expression values for
#' values with each row and column representing genes/isoforms and cells,
#' respectively. Row and column names should be included.
#' Expression values should not be log-transformed.
#' @param relative_expr_thresh Relative expression values below this threshold are considered zero.
#'
#' @return an vector of most abundant relative_expr value corresponding to the RPC 1.
#' @export
estimate.t <- function(
    relative_expr_matrix,
    relative_expr_thresh = 0.1) {
  # apply each column
  unlist(
    apply(
      relative_expr_matrix, 2, function(relative_expr) {
        10^mean(transcript.mode(log10(relative_expr[relative_expr > relative_expr_thresh])))
      }
    )
  ) # avoid multiple output
}

#' use gaussian kernel to calculate the mode of transcript counts
#'
#' @param x log tranformed relative expression
#' @param  breaks control parameter
#'
#' @export
transcript.mode <- function(
    x,
    breaks = "Sturges") {
  if (length(x) < 2) {
    return(0)
  }
  den <- stats::density(x, kernel = c("gaussian"))
  (den$x[den$y == max(den$y)])
}
