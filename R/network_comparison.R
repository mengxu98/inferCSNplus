# =================== Functions to compare networks =====================

#' Function to bootstrap Epoch reconstruction
#'
#' @param matrix genes-by-cells expression matrix
#' @param meta_data sample table
#' @param tfs TFs
#' @param n_reconstructions the number of networks to reconstruct
#' @param ncells_per_sample the number of cells to sample for each reconstruction
#' @param p_value p-value threshold for selecting dynamic genes
#' @param zThresh z-score threshold in static reconstruction
#' @param pseudotime_column column name in meta_data with pseudotime annotation
#' @param group_column column name in meta_data for restricting reconstruction
#' @param method CLR method, either "pearson" or "MI"
#' @param crossweight whether or not to perform crossweighting
#' @param limit_to vector of genes. If not NULL, will skip finding dynamic genes and use genes in limit_to for reconstruction
#'
#' @return list of grnDFs
#'
#' @export
sample_and_epoch_reconstruct <- function(
    matrix,
    meta_data,
    tfs,
    n_reconstructions = 10,
    ncells_per_sample = 400,
    p_value = 0.05,
    zThresh = 5,
    pseudotime_column = "pseudotime",
    group_column = "leiden",
    method = "pearson",
    crossweight = FALSE,
    limit_to = NULL) {
  res <- list()

  for (i in 1:n_reconstructions) {
    print(i)
    # sample
    meta_data_subset <- meta_data[sample(rownames(meta_data), ncells_per_sample, replace = FALSE), ]
    matrix_subset <- matrix[, rownames(meta_data_subset)]
    matrix_subset <- matrix_subset[rowSums(matrix_subset) != 0, ]
    tfs_subset <- intersect(tfs, rownames(matrix_subset))

    # reconstruct
    if (is.null(limit_to)) {
      xdyn <- findDynGenes(matrix_subset, meta_data_subset, unique(meta_data_subset[, group_column]), group_column = group_column, pseudotime_column = pseudotime_column)
      dgenes <- names(xdyn$genes)[xdyn$genes < p_value]
      dgenes <- dgenes[!is.na(dgenes)]
    } else {
      dgenes <- limit_to
      dgenes <- intersect(dgenes, rownames(matrix_subset))
    }

    grnDF <- inferCSN(matrix_subset[dgenes, ], regulators = tfs_subset)

    if (crossweight) {
      meta_data_subset <- meta_data_subset[order(meta_data_subset[, pseudotime_column], decreasing = FALSE), ]
      grnDF <- crossweight(grnDF, matrix_subset[, rownames(meta_data_subset)])
    }

    res[[i]] <- grnDF
  }

  res
}

# =================== Differential Network Functions =====================
#' Function to compute edge differences between networks
#'
#' @param grnDFs list of static grnDFs
#' @param tfs TFs
#' @param weight_column column name in grnDFs containing edge weight
#'
#' @return dataframe of edges with difference scores between the GRNs
#'
#' @export
edge_uniqueness <- function(
    grnDFs,
    tfs,
    weight_column) {
  # all genes
  genes <- unlist(lapply(grnDFs, function(x) union(x$target, x$regulator)), use.names = FALSE)
  genes <- unique(genes)
  tfs <- intersect(genes, tfs)

  # cast grnDFs into adjacency matrices -- have to do this rather than just
  # using dataframes... roundabout, but allows us to keep track of 0-weighted edges too...
  adj_list <- sapply(names(grnDFs), function(x) NULL)

  for (df in names(adj_list)) {
    graph <- igraph::graph_from_data_frame(grnDFs[[df]][, c("regulator", "target", weight_column)], directed = TRUE)
    addvtcs <- setdiff(genes, igraph::V(graph)$name)
    graph <- igraph::add_vertices(graph, length(addvtcs), name = addvtcs)
    adj <- igraph::as_adjacency_matrix(graph, attr = weight_column)

    # prune matrix
    adj <- adj[tfs, ]
    adj <- t(as.matrix(adj))

    adj_list[[df]] <- adj
  }

  # for each edge, compute max-min edge weight (highest scoring treatment/condition - lowest).
  # This should pull out the most variable edges by magnitude.
  # Will need to check how many of these are due to nodes not in network.
  # compare against random?

  score_df <- edge_rank(adj_list)
  score_df
}

#' Function to compute edge differences
#'
#' @param grnMats list of adjacency matrices
#'
#' @return dataframe of edges with difference scores between the GRNs
edge_rank <- function(grnMats) {
  full_df <- as.data.frame(as.table(as.matrix(grnMats[[1]])))
  full_df <- full_df[, 1:2]
  colnames(full_df) <- c("target", "regulator")

  for (df in names(grnMats)) {
    add <- as.data.frame(as.table(as.matrix(grnMats[[df]])))
    colnames(add) <- c("target", "regulator", df)

    full_df <- merge(full_df, add, by = c("target", "regulator"))
  }

  # compute min and max scores; diff score
  full_df$min <- apply(full_df[, !names(full_df) %in% c("target", "regulator")], 1, FUN = min)
  full_df$max <- apply(full_df[, !names(full_df) %in% c("target", "regulator")], 1, FUN = max)

  full_df$diff <- full_df$max - full_df$min

  # order
  full_df <- full_df[order(full_df$diff, decreasing = TRUE), ]

  full_df
}


# edgeDF is the result of running 'edge_uniqueness' followed by split into dynamic network
# type DE_on vs DE_off network == "on" or "off"

#' Compute a dynamic difference network
#'
#'
#'
#' @param edgeDF the result of running edge_uniqueness
#' @param epochs list of epoch gene assignments
#' @param condition condition of interest, should be one of the treatment (aka network name) or column names in edgeDF
#' @param type "on" or "off" specifies either finding edges that are active in the condition network but off others or inactive in condition network but active in others
#' @param diff_thresh edge difference threshold to determine if edge is uniquely on or off
#' @param condition_thresh edge weight theshold in condition network to keep or filter out in difference network
#'
#' @return list of dataframes representing the dynamic difference network
#'
#' @export
#'
dynamic_difference_network <- function(
    edgeDF,
    epochs,
    condition,
    type,
    diff_thresh = 3,
    condition_thresh = 6) {
  edgeDF <- edgeDF[edgeDF$diff != 0, ]

  # sort edges by epoch (limit to dynamic, union of all networks)
  epochs <- lapply(epochs, function(x) {
    x$mean_expression <- NULL
    x
  })

  temp <- vector(mode = "list", length = length(names(epochs[[1]])))
  names(temp) <- names(epochs[[1]])

  for (e in names(epochs)) {
    temp <- Map(c, temp, epochs[[e]])
  }

  epochs <- sapply(temp, unique)

  dynamic_edges <- list()
  for (epoch in names(epochs)) {
    dynamic_edges[[epoch]] <- edgeDF[edgeDF$regulator %in% epochs[[epoch]], ]
  }

  edgeDF <- dynamic_edges


  conditions <- names(edgeDF[[1]])[!(names(edgeDF[[1]]) %in% c("target", "regulator", "min", "max", "diff"))]
  if (!(condition %in% conditions)) {
    stop("condition not in network.")
  }
  diffnet <- lapply(edgeDF, function(x) x[, c("target", "regulator", condition, "diff")])
  diffnet <- lapply(diffnet, function(x) x[x$diff > diff_thresh, ])

  if (type == "on") {
    diffnet <- Map(function(x, y) x[x$regulator %in% y, ], diffnet, epochs)
    diffnet <- lapply(diffnet, function(x) x[x[, condition] >= condition_thresh, ])
  } else if (type == "off") {
    diffnet <- lapply(diffnet, function(x) x[x[, condition] < condition_thresh, ])
  } else {
    stop("type should be 'on' or 'off'.")
  }

  diffnet
}


#' Perform community detection on a dynamic network
#'
#' @param diffnet diffnet
#' @param method community detection method, currently only "louvain"
#' @param use_weights whether or not to use edge weights (for weighted graphs)
#' @param weight_column if using weights, name of the column containing edge weights
#'
#' @return community assignments of nodes in the dynamic network
#'
#' @export
diffnet_community_detection <- function(
    diffnet,
    method = "louvain",
    use_weights = FALSE,
    weight_column = NULL) {
  # if directed=TRUE, remember to flip target and regulator in diffnet dfs
  graphs <- lapply(diffnet, function(x) {
    g <- igraph::graph_from_data_frame(x, directed = FALSE)
    g
  })

  if (use_weights) {
    weights <- igraph::edge_attr(x, weight_column)
  } else {
    weights <- NA
  }

  if (method == "louvain") {
    communities <- lapply(graphs, function(x) {
      c <- igraph::cluster_louvain(x, weights = weights)
      c
    })
  }

  communities
}

# useful function to add correlation-based interaction type to diffnet after using dynamic_difference_network
#' Adds interaction type to dynamic differential network
#'
#' @param diffnet diffnet
#' @param type "on" or "off" depending on the type of differential network. If "on" will assign type based on grnDF_on. Otherwise interaction assigned from grnDF_offlist.
#' @param grnDF_on the static network in which diffnet edges are active
#' @param grnDF_offlist list of static networks in which diffnet edges are inactive
#'
#' @return community assignments of nodes in the dynamic network
#'
#' @export
add_type <- function(
    diffnet,
    type,
    grnDF_on,
    grnDF_offlist) {
  diffnet <- diffnet[sapply(diffnet, function(x) dim(x)[1]) > 0]
  fun <- function(x, y) {
    merge(x, y, by = c("target", "regulator"), all.x = TRUE)
  }
  if (type == "on") {
    added <- lapply(diffnet, fun, grnDF_on[, c("target", "regulator", "corr")])
    added <- lapply(added, transform, interaction = ifelse(corr < 0, "repression", "activation"))
  } else if (type == "off") {
    added <- diffnet
    for (i in 1:length(grnDF_offlist)) {
      added <- lapply(added, fun, grnDF_offlist[[i]][, c("target", "regulator", "corr")])
      # added<-lapply(added,function(x) {colnames(x)[colnames(x)=="corr"]<-paste0("othernet",i); x})
    }
    added <- lapply(added, function(x) {
      cols <- colnames(x)[grepl("corr", colnames(x))]
      x$interaction <- NA
      x[(rowSums((x[, cols] >= 0) | (is.na(x[, cols])), na.rm = TRUE)) == length(cols), "interaction"] <- "activation"
      x[(rowSums((x[, cols] < 0) | (is.na(x[, cols])), na.rm = TRUE)) == length(cols), "interaction"] <- "repression"
      x
    })
  } else {
    stop("invalid type.")
  }

  added
}


# =================== General, Broad comparison functions=====================

#' Computes frobenius distance in a pairwise manner between two sets of networks
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param weight_column	column name containing edge weights
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of frobenius distances
#'
#' @export
compute_frobenius_distance <- function(
    netlist1,
    netlist2,
    weight_column = "weight",
    compare_within_netlist1 = TRUE,
    compare_within_netlist2 = TRUE) {
  # between networks in netlist1
  if (compare_within_netlist1) {
    df1 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()
    for (i in 1:nrow(df1)) {
      net1 <- netlist1[[df1$X1[i]]]
      net2 <- netlist1[[df1$X2[i]]]

      net1 <- reshape2::acast(
        net1,
        regulator ~ target,
        value.var = weight_column
      )
      net1[is.na(net1)] <- 0
      net1 <- t(net1)

      net2 <- reshape2::acast(
        net2,
        regulator ~ target,
        value.var = weight_column
      )
      net2[is.na(net2)] <- 0
      net2 <- t(net2)

      rows <- union(rownames(net1), rownames(net2))
      cols <- union(colnames(net1), colnames(net2))

      # theres a better solution than this....
      missing <- setdiff(rows, rownames(net1))
      addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
      rownames(addnet1) <- missing
      colnames(addnet1) <- colnames(net1)
      net1 <- rbind(net1, addnet1)

      missing <- setdiff(rows, rownames(net2))
      addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
      rownames(addnet2) <- missing
      colnames(addnet2) <- colnames(net2)
      net2 <- rbind(net2, addnet2)

      missing <- setdiff(cols, colnames(net1))
      addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
      rownames(addnet1) <- rownames(net1)
      colnames(addnet1) <- missing
      net1 <- cbind(net1, addnet1)

      missing <- setdiff(cols, colnames(net2))
      addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
      rownames(addnet2) <- rownames(net2)
      colnames(addnet2) <- missing
      net2 <- cbind(net2, addnet2)

      # order rows and columns, take difference, then norm
      net1 <- net1[rows, cols]
      net2 <- net2[rows, cols]

      diff <- net1 - net2
      score <- c(score, norm(diff, "F"))
    }

    df1$score <- score
  }

  # between networks in netlist2
  if (compare_within_netlist2) {
    df2 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()

    for (i in 1:nrow(df2)) {
      net1 <- netlist2[[df2$X1[i]]]
      net2 <- netlist2[[df2$X2[i]]]

      net1 <- reshape2::acast(
        net1,
        regulator ~ target,
        value.var = weight_column
      )
      net1[is.na(net1)] <- 0
      net1 <- t(net1)

      net2 <- reshape2::acast(
        net2,
        regulator ~ target,
        value.var = weight_column
      )
      net2[is.na(net2)] <- 0
      net2 <- t(net2)

      rows <- union(rownames(net1), rownames(net2))
      cols <- union(colnames(net1), colnames(net2))

      # theres a better solution than this....
      missing <- setdiff(rows, rownames(net1))
      addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
      rownames(addnet1) <- missing
      colnames(addnet1) <- colnames(net1)
      net1 <- rbind(net1, addnet1)

      missing <- setdiff(rows, rownames(net2))
      addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
      rownames(addnet2) <- missing
      colnames(addnet2) <- colnames(net2)
      net2 <- rbind(net2, addnet2)

      missing <- setdiff(cols, colnames(net1))
      addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
      rownames(addnet1) <- rownames(net1)
      colnames(addnet1) <- missing
      net1 <- cbind(net1, addnet1)

      missing <- setdiff(cols, colnames(net2))
      addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
      rownames(addnet2) <- rownames(net2)
      colnames(addnet2) <- missing
      net2 <- cbind(net2, addnet2)

      # order rows and columns, take difference, then norm
      net1 <- net1[rows, cols]
      net2 <- net2[rows, cols]

      diff <- net1 - net2
      score <- c(score, norm(diff, "F"))
    }

    df2$score <- score
  }

  # between networks in netlist1 and netlist2
  df3 <- data.frame(
    expand.grid(1:length(netlist1), 1:length(netlist2))
  )
  colnames(df3) <- c("X1", "X2")
  score <- c()

  for (i in 1:nrow(df3)) {
    net1 <- netlist1[[df3$X1[i]]]
    net2 <- netlist2[[df3$X2[i]]]

    net1 <- reshape2::acast(
      net1,
      regulator ~ target,
      value.var = weight_column
    )
    net1[is.na(net1)] <- 0
    net1 <- t(net1)

    net2 <- reshape2::acast(
      net2,
      regulator ~ target,
      value.var = weight_column
    )
    net2[is.na(net2)] <- 0
    net2 <- t(net2)

    rows <- union(rownames(net1), rownames(net2))
    cols <- union(colnames(net1), colnames(net2))

    # theres a better solution than this....
    missing <- setdiff(rows, rownames(net1))
    addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
    rownames(addnet1) <- missing
    colnames(addnet1) <- colnames(net1)
    net1 <- rbind(net1, addnet1)

    missing <- setdiff(rows, rownames(net2))
    addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
    rownames(addnet2) <- missing
    colnames(addnet2) <- colnames(net2)
    net2 <- rbind(net2, addnet2)

    missing <- setdiff(cols, colnames(net1))
    addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
    rownames(addnet1) <- rownames(net1)
    colnames(addnet1) <- missing
    net1 <- cbind(net1, addnet1)

    missing <- setdiff(cols, colnames(net2))
    addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
    rownames(addnet2) <- rownames(net2)
    colnames(addnet2) <- missing
    net2 <- cbind(net2, addnet2)

    # order rows and columns, take difference, then norm
    net1 <- net1[rows, cols]
    net2 <- net2[rows, cols]

    diff <- net1 - net2
    score <- c(score, norm(diff, "F"))
  }

  df3$score <- score

  # summarize
  if (compare_within_netlist1) {
    df1$group <- "netlist1"
  }
  if (compare_within_netlist2) {
    df2$group <- "netlist2"
  }
  df3$group <- "cross_comparison"

  if (compare_within_netlist1 & compare_within_netlist2) {
    df <- rbind(df1, df2)
    df <- rbind(df, df3)
  } else if (compare_within_netlist2) {
    df <- rbind(df2, df3)
  } else if (compare_within_netlist1) {
    df <- rbind(df1, df3)
  } else {
    df <- df3
  }

  df
}

# compute Jaccard index of top regulators
# right now only PageRank supported

#' Computes Jaccard similarity between top regulators in two sets of networks
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param n_regs the number of top regulators to compare from each network
#' @param method method to find top regulators. Currently only supports "pagerank"
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of Jaccard similarities of top regulators
#'
#' @export
#'
compute_JI_topregs <- function(
    netlist1,
    netlist2,
    n_regs = 15,
    method = "pagerank",
    compare_within_netlist1 = TRUE,
    compare_within_netlist2 = TRUE) {
  # between networks in netlist1
  if (compare_within_netlist1) {
    df1 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()

    for (i in 1:nrow(df1)) {
      net1 <- netlist1[[df1$X1[i]]]
      net2 <- netlist1[[df1$X2[i]]]

      net1_ranks <- compute_pagerank(list(X = net1))
      net2_ranks <- compute_pagerank(list(X = net2))

      net1_topregs <- net1_ranks$X$gene[1:n_regs]
      net2_topregs <- net2_ranks$X$gene[1:n_regs]

      ji <- length(
        intersect(net1_topregs, net2_topregs)
      ) / length(union(net1_topregs, net2_topregs))

      score <- c(score, ji)
    }

    df1$jaccard <- score
  }

  # between networks in netlist2
  if (compare_within_netlist2) {
    df2 <- data.frame(t(utils::combn(1:length(netlist2), 2)))
    score <- c()

    for (i in 1:nrow(df2)) {
      net1 <- netlist2[[df2$X1[i]]]
      net2 <- netlist2[[df2$X2[i]]]

      net1_ranks <- compute_pagerank(list(X = net1))
      net2_ranks <- compute_pagerank(list(X = net2))

      net1_topregs <- net1_ranks$X$gene[1:n_regs]
      net2_topregs <- net2_ranks$X$gene[1:n_regs]

      ji <- length(
        intersect(net1_topregs, net2_topregs)
      ) / length(union(net1_topregs, net2_topregs))

      score <- c(score, ji)
    }

    df2$jaccard <- score
  }

  # between networks in netlist1 and netlist2
  df3 <- data.frame(expand.grid(1:length(netlist1), 1:length(netlist2)))
  colnames(df3) <- c("X1", "X2")
  score <- c()

  for (i in 1:nrow(df3)) {
    net1 <- netlist1[[df3$X1[i]]]
    net2 <- netlist2[[df3$X2[i]]]

    net1_ranks <- compute_pagerank(list(X = net1))
    net2_ranks <- compute_pagerank(list(X = net2))

    net1_topregs <- net1_ranks$X$gene[1:n_regs]
    net2_topregs <- net2_ranks$X$gene[1:n_regs]

    ji <- length(
      intersect(net1_topregs, net2_topregs)
    ) / length(union(net1_topregs, net2_topregs))

    score <- c(score, ji)
  }

  df3$jaccard <- score

  # summarize
  if (compare_within_netlist1) {
    df1$group <- "netlist1"
  }
  if (compare_within_netlist2) {
    df2$group <- "netlist2"
  }
  df3$group <- "cross_comparison"

  if (compare_within_netlist1 & compare_within_netlist2) {
    df <- rbind(df1, df2)
    df <- rbind(df, df3)
  } else if (compare_within_netlist2) {
    df <- rbind(df2, df3)
  } else if (compare_within_netlist1) {
    df <- rbind(df1, df3)
  } else {
    df <- df3
  }

  df
}


#' Computes Jaccard similarity between top regulators in two sets of networks across a range of top X regulators
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param n_regs a vector indicating which values of top regulators to scan across
#' @param func func
#' @param method method to find top regulators. Currently only supports "pagerank"
#' @param weight_column column name in grnDFs containing edge weights
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of Jaccard similarities of top regulators
#'
#' @export
JI_across_topregs <- function(
    netlist1,
    netlist2,
    n_regs = 3:15,
    func = "mean",
    method = "pagerank",
    weight_column = "zscore",
    compare_within_netlist1 = TRUE,
    compare_within_netlist2 = TRUE) {
  if (func == "mean") {
    res <- data.frame(
      group = character(),
      jaccard = numeric(),
      n_regs = numeric()
    )
    for (i in n_regs) {
      ji <- compute_JI_topregs(
        netlist1,
        netlist2,
        n_regs = i,
        method = method,
        compare_within_netlist1 = compare_within_netlist1,
        compare_within_netlist2 = compare_within_netlist2
      )
      ji <- ji[, c("group", "jaccard")]
      ji <- aggregate(. ~ group, ji, mean)
      ji$n_regs <- i

      res <- rbind(res, ji)
    }
  }
  res
}

tally_topregs <- function(
    netlist1,
    netlist2,
    n_regs = 15,
    method = "pagerank") {
  # in netlist1
  topregs_1 <- c()
  for (i in 1:length(netlist1)) {
    net1 <- netlist1[[i]]
    net1_ranks <- compute_pagerank(list(X = net1))
    net1_topregs <- net1_ranks$X$gene[1:n_regs]

    topregs_1 <- c(topregs_1, net1_topregs)
  }

  # in netlist2
  topregs_2 <- c()
  for (i in 1:length(netlist2)) {
    net2 <- netlist2[[i]]
    net2_ranks <- compute_pagerank(list(X = net2))
    net2_topregs <- net2_ranks$X$gene[1:n_regs]

    topregs_2 <- c(topregs_2, net2_topregs)
  }

  # summarize
  topregs_1 <- table(topregs_1)
  topregs_2 <- table(topregs_2)

  list(netlist1 = topregs_1, netlist2 = topregs_2)
}

average_pagerank_ranks <- function(
    netlist1,
    max_rank = NULL) {
  if (is.null(names(netlist1))) {
    names(netlist1) <- paste0("X", seq(1:length(netlist1)))
  }

  if (length(netlist1) == 1) {
    net_ranks <- compute_pagerank(netlist1)
    net_ranks <- lapply(net_ranks, function(x) x[x$is_regulator, ])
    net_ranks <- lapply(net_ranks, function(x) cbind(x, rank = 1:nrow(x)))
    net_ranks <- lapply(net_ranks, function(x) x[, c("gene", "rank")])
    res <- net_ranks[[1]]
    colnames(res) <- c("gene", "mean_rank")
    return(res)
  }

  net_ranks <- compute_pagerank(netlist1)
  net_ranks <- lapply(net_ranks, function(x) x[x$is_regulator, ])
  net_ranks <- lapply(net_ranks, function(x) cbind(x, rank = 1:nrow(x)))
  net_ranks <- lapply(net_ranks, function(x) x[, c("gene", "rank")])
  net_ranks <- lapply(net_ranks, function(x) as.data.frame(t(x)))

  if (is.null(max_rank)) {
    max_rank <- max(unlist(lapply(net_ranks, function(x) nrow(x))))
  }

  ranks <- plyr::rbind.fill(net_ranks, fill = NA)

  ranks <- ranks[seq(2, nrow(ranks), 2), ]

  ranks <- data.frame(sapply(ranks, function(x) as.numeric(as.character(x))))
  ranks[is.na(ranks)] <- max_rank
  avgranks <- colMeans(ranks)

  res <- data.frame(gene = names(avgranks), mean_rank = avgranks)
  res <- res[order(res$mean_rank, decreasing = FALSE), ]

  res
}


# function that takes a network as grnDF, returns weighted betweennes of TFs
compute_betweenness <- function(
    grnDF,
    tfs = NULL,
    normalized = TRUE) {
  if (is.null(tfs)) {
    tfs <- unique(grnDF$regulator)
  }

  grnDF <- grnDF[, c("regulator", "target", "weight")]
  g <- igraph::graph_from_data_frame(grnDF, directed = FALSE)

  betweenness <- igraph::betweenness(
    g,
    tfs,
    directed = FALSE,
    normalized = normalized
  )

  betweenness
}

# function that takes a list of network dataframes
# returns long dataframe listing network, betweenness of each regulator, degree of each regulator

#' Computes betweenness and degree of each regulator for each network in a list of networks
#'
#' @param netlist netlist
#' @param weight_column column name in grnDFs containing edge weights
#' @param normalized whether or not to normalize degree and betweenness
#'
#' @return dataframe listing network, betweenness of each regulator, degree of each regulator
#'
#' @export
biglist_compute_betweenness_degree <- function(
    netlist,
    weight_column = "zscore",
    normalized = TRUE) {
  res <- data.frame(
    betweenness = numeric(),
    degree = numeric(),
    network = character(),
    regulator = character()
  )

  if (is.null(names(netlist))) {
    names(netlist) <- seq(1:length(netlist))
  }

  for (name in names(netlist)) {
    net <- netlist[[name]]
    g <- igraph::graph_from_data_frame(
      net[, c("regulator", "target")],
      directed = FALSE
    )
    betweenness <- compute_betweenness(
      net,
      normalized = normalized
    )
    degree <- degree(
      g,
      unique(net$regulator),
      mode = "all",
      normalized = normalized
    )

    df <- cbind(betweenness, degree = degree[names(betweenness)])
    df <- as.data.frame(df)
    df$network <- name
    df$regulator <- rownames(df)

    res <- rbind(res, df)
  }

  res
}


# =================== Useful Plotting Functions =====================

#' Plot the dynamic differential network
#'
#' @param grn the dynamic network
#' @param tfs TFs
#' @param only_TFs whether or not to only plot TFs and exclude non-regulators
#' @param order the order in which to plot epochs, or which epochs to plot
#'
#' @return plot
#'
#' @export
#'
plot_dyn_diffnet <- function(
    grn,
    tfs,
    only_TFs = TRUE,
    order = NULL) {
  g <- list()

  if (!is.null(order)) {
    grn <- grn[order]
  }

  for (i in 1:length(grn)) {
    df <- grn[[i]]

    if (only_TFs) {
      df <- df[df$target %in% tfs, ]
    }
    if (nrow(df) == 0) {
      next
    }

    net <- igraph::graph_from_data_frame(
      df[, c("regulator", "target", "interaction")],
      directed = FALSE
    )
    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    tfnet <- ggnetwork::ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    tfnet$is_regulator <- as.character(tfnet$name %in% tfs)

    cols <- c("activation" = "blue", "repression" = "red")
    g[[i]] <- ggplot() +
      geom_edges(
        data = tfnet,
        aes(x = x, y = y, xend = xend, yend = yend, color = interaction),
        size = 0.75,
        curvature = 0.1,
        alpha = .6
      ) +
      scale_color_manual(values = cols) +
      geom_nodes(
        data = tfnet,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "darkgray",
        size = 6,
        alpha = .5
      ) +
      geom_nodes(
        data = tfnet[tfnet$is_regulator == "TRUE", ],
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "#8C4985",
        size = 6,
        alpha = .8
      ) +
      geom_nodelabel_repel(
        data = tfnet,
        aes(x = x, y = y, label = name),
        size = 2.5,
        color = "#5A8BAD"
      ) +
      theme_blank() +
      ggtitle(names(grn)[i])

    g[[i]] <- g[[i]] + theme(legend.position = "none")
  }
  g <- g[!sapply(g, is.null)]
  do.call(gridExtra::grid.arrange, g)
}


# same as above, but will also do community detection and betweenness computation, and color/fade accordingly

#' Plot the dynamic differential network but colored by communities and optionally faded by betweenness
#'
#' @param grn the dynamic network
#' @param tfs TFs
#' @param only_TFs whether or not to only plot TFs and exclude non-regulators
#' @param order the order in which to plot epochs, or which epochs to plot
#' @param compute_betweenness whether or not to fade nodes by betweenness
#'
#' @return plot
#'
#' @export
plot_diffnet_detail <- function(
    grn,
    tfs,
    only_TFs = TRUE,
    order = NULL,
    compute_betweenness = TRUE) {
  g <- list()

  if (!is.null(order)) {
    grn <- grn[order]
  }

  if (any(sapply(grn, function(x) {
    ("interaction" %in% names(x))
  }) == FALSE)) {
    grn <- add_interaction_type(grn)
  }

  for (i in 1:length(grn)) {
    df <- grn[[i]]

    if (only_TFs) {
      df <- df[df$target %in% tfs, ]
    }
    if (nrow(df) == 0) {
      next
    }

    df <- df[, c("regulator", "target", "weight", "interaction")]

    net <- igraph::graph_from_data_frame(df, directed = FALSE)

    if (compute_betweenness) {
      b <- igraph::betweenness(net, directed = FALSE, normalized = TRUE)
      b <- as.data.frame(b)
      colnames(b) <- "betweenness"
      b$gene <- rownames(b)
      b <- b[, c("gene", "betweenness")]
      c <- as.data.frame(as.table(igraph::membership(igraph::cluster_louvain(net))))
      colnames(c) <- c("gene", "communities")
      vtx_features <- merge(b, c, by = "gene", all = TRUE)
    } else {
      c <- as.data.frame(as.table(igraph::membership(igraph::cluster_louvain(net))))
      colnames(c) <- c("gene", "communities")
      vtx_features <- c
    }

    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    tfnet <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    tfnet$is_regulator <- as.character(tfnet$name %in% tfs)
    if (compute_betweenness) {
      tfnet$betweenness <- vtx_features$betweenness[match(tfnet$name, vtx_features$gene)]
    }
    tfnet$communities <- as.factor(vtx_features$communities[match(tfnet$name, vtx_features$gene)])

    cols <- c("activation" = "blue", "repression" = "red")
    num_cols2 <- length(unique(tfnet$communities))
    if (num_cols2 <= 8) {
      cols2 <- RColorBrewer::brewer.pal(num_cols2, "Set1")
    } else {
      cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_cols2)
    }
    names(cols2) <- unique(tfnet$communities)
    cols <- c(cols, cols2)

    if (compute_betweenness) {
      g[[i]] <- ggplot() +
        geom_edges(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend, color = interaction), size = 0.75, curvature = 0.1, alpha = .6) +
        scale_color_manual(values = cols) +
        geom_nodes(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgray", size = 6, alpha = .5) +
        geom_nodes(data = tfnet[tfnet$is_regulator == "TRUE", ], aes(x = x, y = y, xend = xend, yend = yend, color = communities, alpha = betweenness + .5), size = 6) +
        geom_nodelabel_repel(data = tfnet, aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD") +
        theme_blank() +
        ggtitle(names(grn)[i])
    } else {
      g[[i]] <- ggplot() +
        geom_edges(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend, color = interaction), size = 0.75, curvature = 0.1, alpha = .6) +
        scale_color_manual(values = cols) +
        geom_nodes(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgray", size = 6, alpha = .5) +
        geom_nodes(data = tfnet[tfnet$is_regulator == "TRUE", ], aes(x = x, y = y, xend = xend, yend = yend, color = communities), alpha = .5, size = 6) +
        geom_nodelabel_repel(data = tfnet, aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD") +
        theme_blank() +
        ggtitle(names(grn)[i])
    }

    g[[i]] <- g[[i]] + theme(legend.position = "none")
  }
  g <- g[!sapply(g, is.null)]
  do.call(gridExtra::grid.arrange, g)
}
