#' @title Detection of metacells from single-cell gene expression matrix
#'
#' @param matrix log-normalized gene expression matrix with rows to be genes and cols to be cells
#' @param genes_use a vector of genes used to compute PCA
#' @param genes_exclude a vector of genes to be excluded when computing PCA
#' @param var_genes_num if \code{"genes.use"} is not provided, \code{"var_genes_num"} genes with the largest variation are used
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)
#' @param knn_k parameter to compute single-cell kNN network
#' @param do_scale whether to scale gene expression matrix when computing PCA
#' @param pc_num number of principal components to use for construction of single-cell kNN network
#' @param fast_pca use \link[irlba]{irlba} as a faster version of prcomp
#' @param do_approx compute approximate kNN in case of a large dataset (>50'000)
#' @param approx_num number of cells to subsample for an approximate approach
#' @param directed directed
#' @param use_nn2 use_nn2
#' @param seed seed to use to subsample cells for an approximate approach
#' @param weights vector of a cell weight (NULL by default), used for computing average gene expression withing cluster of metaells
#' @param do_median_norm whether to normalize by median value (FALSE by default)
#' @param block_size number of cells to map to the nearest metacell at the time (for approx coarse-graining)
#' @param cluster_method clustering method to identify metacells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).
#' @param ... Parameters for other methods.
#'
#' @return a metacell matrix
#'
#' @references
#' https://github.com/GfellerLab/SuperCell
#'
#' https://github.com/kuijjerlab/SCORPION
#'
#' @export
#'
#' @examples
#' data("example_matrix")
#' meta_cells_matrix <- meta_cells(
#'   example_matrix,
#'   pc_num = 10,
#'   fast_pca = FALSE
#' )
#' dim(meta_cells_matrix)
#' meta_cells_matrix[1:6, 1:6]
meta_cells <- function(
    matrix,
    genes_use = NULL,
    genes_exclude = NULL,
    var_genes_num = min(1000, nrow(matrix)),
    gamma = 10,
    knn_k = 5,
    do_scale = TRUE,
    pc_num = 25,
    fast_pca = TRUE,
    do_approx = FALSE,
    approx_num = 20000,
    directed = FALSE,
    use_nn2 = TRUE,
    seed = 1,
    cluster_method = "walktrap",
    block_size = 10000,
    weights = NULL,
    do_median_norm = FALSE,
    ...) {
  matrix <- t(matrix)
  cells_num <- ncol(matrix)
  matrix_raw <- matrix
  cluster_method <- match.arg(
    cluster_method,
    c("walktrap", "louvain")
  )

  matrix <- matrix[rowSums(matrix) > 0, ]
  keep.genes <- setdiff(rownames(matrix), genes_exclude)
  matrix <- matrix[keep.genes, ]

  if (is.null(genes_use)) {
    var_genes_num <- min(var_genes_num, nrow(matrix))
    if (cells_num > 50000) {
      set.seed(seed)
      idx <- sample(cells_num, 50000)
      gene_var <- apply(matrix[, idx], 1, stats::var)
    } else {
      gene_var <- apply(matrix, 1, stats::var)
    }

    genes_use <- names(sort(gene_var, decreasing = TRUE))[1:var_genes_num]
  }

  if (length(intersect(genes_use, genes_exclude)) > 0) {
    stop("Sets of genes_use and genes_exclude have non-empty intersection")
  }

  genes_use <- genes_use[genes_use %in% rownames(matrix)]
  matrix <- matrix[genes_use, ]

  if (do_approx & approx_num >= cells_num) {
    do_approx <- FALSE
    warning(
      "number of obtained metacells is larger or equal to the number of single cells,
    thus, an exact simplification will be performed"
    )
  }

  if (do_approx & (approx_num < round(cells_num / gamma))) {
    approx_num <- round(cells_num / gamma)
    warning("number of obtained metacells is set to ", approx_num)
  }

  if (do_approx & ((cells_num / gamma) > (approx_num / 3))) {
    warning(
      "number of obtained metacells is not much larger than desired number of super-cells,
    so an approximate simplification may take londer than an exact one!"
    )
  }

  if (do_approx) {
    set.seed(seed)
    approx_num <- min(approx_num, cells_num)
    presample <- sample(1:cells_num, size = approx_num, replace = FALSE)
    presampled_cells <- colnames(matrix)[sort(presample)]
    rest_cells <- setdiff(colnames(matrix), presampled_cells)
  } else {
    presampled_cells <- colnames(matrix)
    rest_cells <- c()
  }

  matrix_pca <- Matrix::t(matrix[genes_use, presampled_cells])
  if (do_scale) {
    matrix_pca <- scale(matrix_pca)
  }
  matrix_pca[is.na(matrix_pca)] <- 0

  if (length(pc_num) == 1) {
    pc_num <- 1:pc_num
  }

  if (fast_pca & (cells_num < 1000)) {
    fast_pca <- FALSE
  }

  if (!fast_pca) {
    pca_results <- stats::prcomp(
      matrix_pca,
      rank. = max(pc_num),
      scale. = F,
      center = F
    )
  } else {
    pca_results <- irlba::irlba(matrix_pca, pc_num)
    pca_results$x <- pca_results$u %*% diag(pca_results$d)
    pca_results$rotation <- pca_results$v
  }

  sc_nw <- .build_knn(
    matrix = pca_results$x[, pc_num],
    k = knn_k,
    from = "coordinates",
    use_nn2 = use_nn2,
    dist_method = "euclidean",
    directed = directed
  )

  k <- round(cells_num / gamma)

  if (cluster_method[1] == "walktrap") {
    membership_results <- igraph::cluster_walktrap(
      sc_nw$graph_knn
    ) |>
      igraph::cut_at(k)
  } else if (cluster_method[1] == "louvain") {
    log_message(
      "using ", cluster_method, " method to cluster, gamma is ignored.",
      message_type = "warning"
    )
    membership_results <- igraph::cluster_louvain(
      sc_nw$graph_knn
    )$membership
  }

  names(membership_results) <- presampled_cells

  # SC.NW <- igraph::contract(
  #   sc_nw$graph_knn,
  #   membership_results
  # ) |> igraph::simplify(
  #   remove.loops = T,
  #   edge.attr.comb = "sum"
  # )

  if (do_approx) {
    PCA.averaged.SC <- as.matrix(
      Matrix::t(
        .supercell_ge(
          t(
            pca_results$x[, pc_num]
          ),
          groups = membership_results
        )
      )
    )
    matrix_roration <- Matrix::t(matrix[genes_use, rest_cells])

    if (do_scale) {
      matrix_roration <- scale(matrix_roration)
    }
    matrix_roration[is.na(matrix_roration)] <- 0

    membership_omitted <- c()
    if (is.null(block_size) | is.na(block_size)) {
      block_size <- 10000
    }

    N.blocks <- length(rest_cells) %/% block_size
    if (length(rest_cells) %% block_size > 0) {
      N.blocks <- N.blocks + 1
    }

    if (N.blocks > 0) {
      for (i in 1:N.blocks) {
        # compute knn by blocks
        idx_begin <- (i - 1) * block_size + 1
        idx_end <- min(i * block_size, length(rest_cells))

        cur_rest_cell_ids <- rest_cells[idx_begin:idx_end]

        PCA.ommited <- matrix_roration[cur_rest_cell_ids, ] %*% pca_results$rotation[, pc_num]

        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC)

        membership_omitted_cur <- apply(D.omitted.subsampled, 1, which.min)
        names(membership_omitted_cur) <- cur_rest_cell_ids

        membership_omitted <- c(membership_omitted, membership_omitted_cur)
      }
    }

    membership_all <- c(membership_results, membership_omitted)
    membership_all <- membership_all[colnames(matrix)]
  } else {
    membership_all <- membership_results[colnames(matrix)]
  }

  membership <- membership_all
  matrix <- matrix_raw

  meta_cells_num <- base::as.vector(table(membership))
  j <- rep(1:max(membership), meta_cells_num)

  goups_idx <- base::split(seq_len(ncol(matrix)), membership)
  i <- unlist(goups_idx)

  if (is.null(weights)) {
    matrix_metacells <- matrix %*% Matrix::sparseMatrix(i = i, j = j)
    matrix_metacells <- sweep(matrix_metacells, 2, meta_cells_num, "/")
  } else {
    if (length(weights) != length(membership)) {
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }
    matrix_metacells <- matrix_metacells %*% Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    weighted_supercell_size <- unlist(
      lapply(
        goups_idx,
        FUN = function(x) {
          sum(weights[x])
        }
      )
    )
    matrix_metacells <- sweep(
      matrix_metacells, 2, weighted_supercell_size, "/"
    )
  }

  if (do_median_norm) {
    matrix_metacells <- (matrix_metacells + 0.01) / apply(matrix_metacells + 0.01, 1, stats::median)
  }

  return(
    t(matrix_metacells)
  )
}

.supercell_ge <- function(
    ge,
    groups,
    mode = c("average", "sum"),
    weights = NULL,
    do.median.norm = FALSE) {
  if (ncol(ge) != length(groups)) {
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }

  mode <- mode[1]
  if (!(mode %in% c("average", "sum"))) {
    stop(paste("mode", mode, "is unknown. Available values are 'average' and 'sum'."))
  }

  supercell_size <- as.vector(table(groups))
  j <- rep(1:max(groups), supercell_size)

  goups_idx <- plyr::split_indices(groups)
  i <- unlist(goups_idx)

  if (is.null(weights)) {
    GE.SC <- ge %*% Matrix::sparseMatrix(i = i, j = j)

    if (mode == "average") {
      GE.SC <- sweep(GE.SC, 2, supercell_size, "/")
    }
  } else {
    if (length(weights) != length(groups)) {
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }

    if (mode != "average") {
      stop(paste("weighted averaging is supposted only for mode = 'average', not for", mode))
    }

    GE.SC <- ge %*% Matrix::sparseMatrix(i = i, j = j, x = weights[i])

    weighted_supercell_size <- unlist(
      lapply(
        goups_idx,
        FUN = function(x) {
          sum(weights[x])
        }
      )
    )
    GE.SC <- sweep(GE.SC, 2, weighted_supercell_size, "/")
  }

  if (do.median.norm) {
    GE.SC <- (GE.SC + 0.01) / apply(GE.SC + 0.01, 1, stats::median)
  }

  return(GE.SC)
}

.build_knn <- function(
    matrix,
    k = 5,
    from = c("dist", "coordinates"),
    use_nn2 = TRUE,
    return_neighbors_order = F,
    dist_method = "euclidean",
    cor_method = "pearson",
    p = 2,
    directed = FALSE) {
  method <- match.arg(from, c("dist", "coordinates"))

  if (method == "coordinates") {
    # from coordinates
    if (use_nn2) {
      if (dist_method != "euclidean") {
        stop(
          paste0(
            "Fast nn2 function from RANN package is used, so",
            dist_method,
            "distnce is not acceptable.
        To use nn2 method, please, choose eucleadian distance.
        If you want to use",
            dist_method,
            "distance, please set parameter use_nn2 to FALSE"
          )
        )
      }
      mode <- ifelse(directed, "out", "all")
      return(
        .build_nn2(matrix = matrix, k = k, mode = mode)
      )
    } else {
      dist_method_ <- match.arg(
        dist_method,
        c(
          "cor",
          "euclidean",
          "maximum",
          "manhattan",
          "canberra",
          "binary",
          "minkowski"
        )
      )

      if (dist_method_ == "cor") {
        cor_method_ <- match.arg(
          cor_method,
          c("pearson", "kendall", "spearman")
        )

        matrix <- stats::as.dist(
          as.matrix(
            1 - stats::cor(t(matrix), method = cor_method)
          )
        )
      } else {
        matrix <- stats::dist(matrix, method = dist_method)
      }
    }
  } else {
    if (use_nn2) {
      stop(
        "Method nn2 cannot be applied to distance, to use fast nn2 method,
      please provide coordinates rather than distance and set parameter from to coordinates"
      )
    }
    return(
      .build_knnd(
        D = matrix,
        k = k,
        return_neighbors_order = return_neighbors_order
      )
    )
  }

  # now matrix is distance in any case
  return(
    .build_knnd(
      D = matrix,
      k = k,
      return_neighbors_order = return_neighbors_order
    )
  )
}

.build_knnd <- function(
    D,
    k = 5,
    return_neighbors_order = TRUE,
    mode = "all") {
  if (!methods::is(D, "matrix") | !methods::is(D, "dist")) {
    stop("D (matrix) mast be a matrix or dist!")
  }

  if (!methods::is(D, "dist")) {
    D <- stats::as.dist(D)
  }

  N <- (1 + sqrt(1 + 8 * length(D))) / 2 # number of cells

  if (k >= N) {
    stop("Not enought neighbors in data set!")
  }
  if (k < 1) {
    stop("Invalid number of nearest neighbors, k must be >= 1!")
  }

  row <- function(i, N) {
    return(
      c(
        if (i > 1) {
          D[(i - 1) + c(0:(i - 2)) * (N - 1 - c(1:(i - 1)) / 2)]
        },
        NA,
        if (i < N) {
          D[((i - 1) * (N - 1) - ((i - 1) * (i - 2) / 2) + 1):(((i - 1) * (N - 1) - ((i - 1) * (i - 2) / 2) + 1) + N - i - 1)]
        }
      )
    )
  }

  neighbors <- t(
    sapply(1:N, function(i) {
      order(row(i, N))[1:k]
    })
  )

  adj.knn <- split(
    neighbors,
    rep(
      1:nrow(neighbors),
      times = ncol(neighbors)
    )
  )

  graph_knn <- igraph::graph_from_adj_list(
    adj.knn,
    duplicate = FALSE,
    mode = mode
  )
  graph_knn <- igraph::simplify(
    graph_knn,
    remove.multiple = TRUE
  )
  igraph::E(graph_knn)$weight <- 1

  if (return_neighbors_order) {
    res <- list(
      graph_knn = graph_knn,
      order = neighbors
    )
  } else {
    res <- list(graph_knn = graph_knn)
  }

  return(res)
}

.build_nn2 <- function(
    matrix,
    k = min(5, ncol(matrix)),
    mode = "all") {
  nn2.res <- RANN::nn2(data = matrix, k = k)
  nn2.res <- nn2.res$nn.idx

  graph_knn <- split(
    nn2.res,
    rep(1:nrow(nn2.res), times = ncol(nn2.res))
  ) |>
    igraph::graph_from_adj_list(
      duplicate = FALSE,
      mode = mode
    ) |>
    igraph::simplify(
      remove.multiple = TRUE
    )
  igraph::E(graph_knn)$weight <- 1

  return(list(graph_knn = graph_knn))
}
