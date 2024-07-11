#' @importFrom stats cor median sd var
#' @import Matrix
makeSuperCells <- function(
    X,
    genes.use = NULL,
    genes.exclude = NULL,
    n.var.genes = min(1000, nrow(X)),
    gamma = 10,
    k.knn = 5,
    do.scale = TRUE,
    n.pc = 25,
    fast.pca = TRUE,
    do.approx = FALSE,
    approx.N = 20000,
    directed = FALSE,
    use.nn2 = TRUE,
    seed = 12345,
    igraph.clustering = c("walktrap", "louvain"),
    block.size = 10000,
    weights = NULL,
    do.median.norm = FALSE) {
  X <- t(X)
  N.c <- ncol(X)
  GE <- X

  if (is.null(rownames(X))) {
    if (!(is.null(genes.use) | is.null(genes.exclude))) {
      stop(
        "rownames(X) is Null \nGene expression matrix X is expected to have genes as rownames"
      )
    } else {
      warning(
        "colnames(X) is Null, \nGene expression matrix X is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' "
      )
      rownames(X) <- paste("gene", 1:nrow(X), sep = "_")
    }
  }

  if (is.null(colnames(X))) {
    warning(
      "colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' "
    )
    colnames(X) <- paste("cell", 1:N.c, sep = "_")
  }

  X <- X[rowSums(X) > 0, ]
  keep.genes <- setdiff(rownames(X), genes.exclude)
  X <- X[keep.genes, ]


  if (is.null(genes.use)) {
    n.var.genes <- min(n.var.genes, nrow(X))
    if (N.c > 50000) {
      set.seed(seed)
      idx <- sample(N.c, 50000)
      gene.var <- apply(X[, idx], 1, var)
    } else {
      gene.var <- apply(X, 1, var)
    }

    genes.use <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }

  if (length(intersect(genes.use, genes.exclude)) > 0) {
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }

  genes.use <- genes.use[genes.use %in% rownames(X)]
  X <- X[genes.use, ]

  if (do.approx & approx.N >= N.c) {
    do.approx <- FALSE
    warning(
      "N.approx is larger or equal to the number of single cells, thus, an exact simplification will be performed"
    )
  }

  if (do.approx & (approx.N < round(N.c / gamma))) {
    approx.N <- round(N.c / gamma)
    warning(paste("N.approx is set to N.SC", approx.N))
  }

  if (do.approx & ((N.c / gamma) > (approx.N / 3))) {
    warning(
      "N.approx is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!"
    )
  }

  if (do.approx) {
    set.seed(seed)
    approx.N <- min(approx.N, N.c)
    presample <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids <- c()
  }

  X.for.pca <-
    Matrix::t(X[genes.use, presampled.cell.ids])
  if (do.scale) {
    X.for.pca <- scale(X.for.pca)
  }
  X.for.pca[is.na(X.for.pca)] <- 0

  if (is.null(n.pc[1]) |
      min(n.pc) < 1) {
    stop("Please, provide a range or a number of components to use: n.pc")
  }
  if (length(n.pc) == 1) {
    n.pc <- 1:n.pc
  }

  if (fast.pca & (N.c < 1000)) {
    # warning("Normal PCA is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }

  if (!fast.pca) {
    PCA.presampled <- stats::prcomp(
      X.for.pca,
      rank. = max(n.pc),
      scale. = F,
      center = F
    )
  } else {
    PCA.presampled <- irlba::irlba(X.for.pca, n.pc)
    PCA.presampled$x <- PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }

  sc.nw <- buildKNN(
    X = PCA.presampled$x[, n.pc],
    k = k.knn,
    from = "coordinates",
    use.nn2 = use.nn2,
    dist_method = "euclidean",
    directed = directed
  )

  # simplify

  k <- round(N.c / gamma)

  if (igraph.clustering[1] == "walktrap") {
    g.s <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership <- igraph::cut_at(g.s, k)
  } else if (igraph.clustering[1] == "louvain") {
    warning(paste(
      "igraph.clustering =",
      igraph.clustering,
      ", gamma is ignored"
    ))
    g.s <- igraph::cluster_louvain(sc.nw$graph.knn)
  } else {
    stop(
      paste(
        "Unknown clustering method (",
        igraph.clustering,
        "), please use louvain or walkrtap"
      )
    )
  }

  membership.presampled <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids

  SC.NW <-
    igraph::contract(sc.nw$graph.knn, membership.presampled)
  SC.NW <-
    igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb = "sum")

  if (do.approx) {
    PCA.averaged.SC <- as.matrix(Matrix::t(supercell_GE(t(
      PCA.presampled$x[, n.pc]
    ), groups = membership.presampled)))
    X.for.roration <- Matrix::t(X[genes.use, rest.cell.ids])

    if (do.scale) {
      X.for.roration <- scale(X.for.roration)
    }
    X.for.roration[is.na(X.for.roration)] <- 0

    membership.omitted <- c()
    if (is.null(block.size) | is.na(block.size))
      block.size <- 10000

    N.blocks <- length(rest.cell.ids) %/% block.size
    if (length(rest.cell.ids) %% block.size > 0)
      N.blocks <- N.blocks + 1

    if (N.blocks > 0) {
      for (i in 1:N.blocks) {
        # compute knn by blocks
        idx.begin <- (i - 1) * block.size + 1
        idx.end   <- min(i * block.size,  length(rest.cell.ids))

        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]

        PCA.ommited          <-
          X.for.roration[cur.rest.cell.ids, ] %*% PCA.presampled$rotation[, n.pc]

        D.omitted.subsampled <-
          proxy::dist(PCA.ommited, PCA.averaged.SC)

        membership.omitted.cur        <-
          apply(D.omitted.subsampled, 1, which.min)
        names(membership.omitted.cur) <- cur.rest.cell.ids

        membership.omitted   <-
          c(membership.omitted, membership.omitted.cur)
      }
    }

    membership.all <- c(membership.presampled, membership.omitted)
    membership.all <- membership.all[colnames(X)]
  } else {
    membership.all <- membership.presampled[colnames(X)]
  }

  membership <- membership.all
  X <- GE

  N.SC <- max(membership)
  supercell_size <- base::as.vector(table(membership))
  j <-
    rep(1:N.SC, supercell_size) # column indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  goups.idx <- base::split(seq_len(ncol(X)), membership) # plyr::split_indices(membership)
  i <- unlist(goups.idx) # row indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  if (is.null(weights)) {
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE <- X %*% M.AV
    GE <- sweep(GE, 2, supercell_size, "/")
  } else {
    if (length(weights) != length(membership)) {
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }
    M.AV <- Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    GE <- GE %*% M.AV

    weighted_supercell_size <-
      unlist(lapply(
        goups.idx,
        FUN = function(x) {
          sum(weights[x])
        }
      ))
    GE <- sweep(GE, 2, weighted_supercell_size, "/")
  }

  if (do.median.norm) {
    GE <- (GE + 0.01) / apply(GE + 0.01, 1, median)
  }

  return(t(GE))
}

buildKNN <- function(
    X,
    k = 5,
    from = c("dist", "coordinates"),
    use.nn2 = TRUE,
    return_neighbors_order = F,
    dist_method = "euclidean",
    cor_method = "pearson",
    p = 2,
    directed = FALSE) {
  av.methods <- c("dist", "coordinates")
  method <- pmatch(from[1], av.methods)
  if (is.na(method)) {
    stop(paste(
      "Unlnown method",
      from,
      "Available methods are",
      paste(av.methods, collapse = ", ")
    ))
  }


  if (method == 2) {
    # from coordinates
    if (use.nn2) {
      if (dist_method != "euclidean") {
        stop(
          paste0(
            "Fast nn2 function from RANN package is used, so",
            dist_method,
            "distnce is not acceptable.
            To use nn2 method, please, choose eucleadian distance.
            If you want to use",
            dist_method,
            "distance, please set parameter use.nn2 to FALSE"
          )
        )
      }
      mode <- ifelse(directed, "out", "all")
      return(buildNN2(X = X, k = k, mode = mode))
    } else {
      av.dist <-
        c(
          "cor",
          "euclidean",
          "maximum",
          "manhattan",
          "canberra",
          "binary",
          "minkowski"
        )
      dist_method_ <- pmatch(dist_method, av.dist)
      if (is.na(dist_method_)) {
        stop(
          paste(
            "Unknown distance method:",
            dist_method,
            "Available dist methods are",
            paste(av.dist, collapse = ", ")
          )
        )
      }
      if (dist_method_ == 1) {
        av.cor_methods <- c("pearson", "kendall", "spearman")
        cor_method_ <- pmatch(cor_method, av.cor_methods)
        if (is.na(cor_method_)) {
          stop(
            paste(
              "Unknown cor method:",
              cor_method,
              "Available cor methods are",
              paste(av.cor_methods)
            )
          )
        }
        X <- stats::as.dist(as.matrix(1 - stats::cor(t(X), method = cor_method)))
      } else {
        X <- stats::dist(X, method = dist_method)
      }
    }
  } else {
    if (use.nn2 == TRUE) {
      stop(
        "Method nn2 cannot be applied to distance, to use fast nn2 method, please provide coordinates rather than distance
        and set parameter from to coordinates"
      )
    }
    return(buildKNND(
      D = X,
      k = k,
      return_neighbors_order = return_neighbors_order
    ))
  }

  # now X is distance in any case
  return(buildKNND(
    D = X,
    k = k,
    return_neighbors_order = return_neighbors_order
  ))
}

#' @importFrom methods is
buildKNND <- function(
    D,
    k = 5,
    return_neighbors_order = TRUE,
    mode = "all") {
  if (!is(D, "matrix") | !is(D, "dist")) {
    stop("D (X) mast be a matrix or dist!")
  }

  if (!is(D, "dist")) {
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
    return(c(
      if (i > 1) D[(i - 1) + c(0:(i - 2)) * (N - 1 - c(1:(i - 1)) / 2)],
      NA,
      if (i < N) D[((i - 1) * (N - 1) - ((i - 1) * (i - 2) / 2) + 1):(((i - 1) * (N - 1) - ((i - 1) * (i - 2) / 2) + 1) + N - i - 1)]
    ))
  }

  neighbors <- t(sapply(1:N, function(i) {
    order(row(i, N))[1:k]
  }))

  adj.knn <- split(neighbors, rep(1:nrow(neighbors), times = ncol(neighbors)))

  graph.knn <- igraph::graph_from_adj_list(adj.knn, duplicate = F, mode = mode)
  graph.knn <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

  if (return_neighbors_order) {
    res <- list(
      graph.knn = graph.knn,
      order = neighbors
    )
  } else {
    res <- list(graph.knn = graph.knn)
  }

  return(res)
}

buildNN2 <- function(
    X,
    k = min(5, ncol(X)),
    mode = "all") {
  nn2.res <- RANN::nn2(data = X, k = k)
  nn2.res <- nn2.res$nn.idx

  adj.knn <- split(nn2.res, rep(1:nrow(nn2.res), times = ncol(nn2.res))) # get adj list

  graph.knn <- igraph::graph_from_adj_list(adj.knn, duplicate = F, mode = mode)

  graph.knn <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

  return(list(graph.knn = graph.knn))
}
