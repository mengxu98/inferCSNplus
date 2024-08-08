# =================== Functions to analyze networks =====================

# Updated scoring and target selection
# score_targets adds a column "mean_score" which is just copied from chip-atlas computed average column
# 				adds a column returning the maximum score for a target across samples
# 				adds a column returning the percent frequency of a hit amongst samples based on >threshold

#' Function to score targets of effectors
#'
#' Adds columns for mean binding score, maximum binding score, and percent frequency target is a hit based on threshold
#'
#' @param aList list of dataframes containing binding data for each effector
#' @param threshold binding score threshold to compute hit frequency
#'
#' @return updated list of dataframes
#'
#' @export
score_targets <- function(aList, threshold = 50) {
  # remove STRING score
  aList <- lapply(aList, function(x) {
    x$STRING <- NULL
    x
  })

  # set rownames if not already set
  aList <- lapply(aList, function(x) {
    if (!is.null(x$Target_genes)) {
      rownames(x) <- x$Target_genes
      x$Target_genes <- NULL
    }
    x
  })

  # add mean_score column
  aList <- lapply(aList, function(x) {
    x$mean_score <- x[, grepl("Average", colnames(x))]
    x
  })
  # remove 'average' column
  aList <- lapply(aList, function(x) {
    x <- x[, !grepl("Average", colnames(x))]
  })

  # add max_score
  aList <- lapply(aList, function(x) {
    x$max_score <- apply(x, 1, max)
    x
  })

  # add percent_freq
  aList <- lapply(aList, function(x) {
    columns <- colnames(x)[!(colnames(x) %in% c("Target_genes", "STRING", "mean_score", "max_score"))]
    n_pos <- apply(x[, columns], 1, function(y) {
      sum(y > threshold)
    })
    x$percent_freq <- n_pos / length(columns)
    x
  })

  aList
}

# aList the result of running score_targets
# by_rank if TRUE will return top n_targets targets; if FALSE will use thresholds
# column column name used in either ranking or thresholding the targets
# n_targets the number of targets to return for each effector if by_rank=TRUE

#' Finds binding targets given list of dataframes containing binding info for effectors
#'
#' @param aList the result of running score_targets
#' @param column column name used in either ranking or thresholding the possible targets
#' @param by_rank if TRUE will return top n_targets; if FALSE will use thresholds
#' @param n_targets the number of targets to return for each effector if by_rank=TRUE
#' @param threshold threshold to use if by_rank=FALSE
#'
#' @return list of binding targets for list of effectors
#'
#' @export
find_targets <- function(
    aList,
    column = "max_score",
    by_rank = FALSE,
    n_targets = 2000,
    threshold = NULL) {
  # note: binding scores are -log10(q-value/FDR) from MACS2

  if (by_rank) {
    aList <- lapply(aList, function(x) {
      x <- x[order(x[, column], decreasing = TRUE), ]
      n_targets <- min(n_targets, nrow(x))
      x[1:n_targets, ]
    })
  } else if (!is.null(threshold)) {
    aList <- lapply(aList, function(x) {
      x[(x[, column] > threshold), ]
    })
  } else {
    stop("Insufficient parameters.")
  }

  res <- lapply(aList, function(x) {
    rownames(x)
  })

  res
}


# =================== Functions to find and quantify average subnetwork or module expression across time =====================

# Function to look at average module (group of genes) expression

#' Computes mean expression of groups of genes
#'
#' @param matrix genes-by-cells expression matrix
#' @param module_list list containing grouped genes, each element is a module of genes
#'
#' @return module expression
#'
#' @export
mean_module_expression <- function(
    matrix,
    module_list) {
  res <- data.frame(matrix(ncol = ncol(matrix), nrow = length(module_list)))
  colnames(res) <- colnames(matrix)
  rownames(res) <- names(module_list)
  for (m in names(module_list)) {
    genes <- intersect(module_list[[m]], rownames(matrix))
    exp_sub <- matrix[genes, ]
    if (length(genes) == 1) {
      res[m, ] <- matrix[genes, ]
    } else {
      res[m, ] <- colMeans(exp_sub)
    }
  }
  res
}

#' Computes mean expression of groups of genes in a dynamic network
#'
#' @param matrix genes-by-cells expression matrix
#' @param community_list list of dataframes with community assginemnts. Each dataframe is the result of running subnets.
#'
#' @return subnetwork expression
#'
#' @export
mean_subnetwork_expression <- function(
    matrix,
    community_list) {
  # get list of subnetwork names
  subnets <- c()
  for (e in names(community_list)) {
    subnets <- c(subnets, paste(e, unique(as.character(community_list[[e]]$communities)), sep = "_"))
  }

  # compute subnetwork expression
  res <- data.frame(matrix(ncol = ncol(matrix), nrow = length(subnets)))
  colnames(res) <- colnames(matrix)
  rownames(res) <- subnets
  for (n in subnets) {
    e <- strsplit(n, "_")[[1]][1]
    c <- strsplit(n, "_")[[1]][2]

    genes <- community_list[[e]]$gene[as.character(community_list[[e]]$communities) == as.character(c)]
    exp_sub <- matrix[as.character(genes), ]

    if (length(genes) == 1) {
      res[n, ] <- matrix[as.character(genes), ]
    } else {
      res[n, ] <- colMeans(exp_sub)
    }
  }

  res
}

#' Function to assign nodes to communities via Louvain clustering
#'
#' @param df dataframe containing a static network
#' @param tfs TFs
#' @param tfonly if TRUE will limit network to TFs only
#'
#' @return dataframe containing community assignments for each gene
#'
#' @export
subnets <- function(
    df,
    tfs,
    tfonly = TRUE) {
  if (tfonly) {
    df <- df[df$target %in% tfs, ]
  }

  df <- df[, c("regulator", "target", "weight", "interaction")]

  net <- igraph::graph_from_data_frame(df, directed = FALSE)

  c <- as.data.frame(
    as.table(igraph::membership(igraph::cluster_louvain(net)))
  )
  colnames(c) <- c("gene", "communities")

  c
}

# =================== Functions to analyze flattened networks =====================

# Flatten networks
# function for within an epoch network, or within a static network
flatten_network <- function(
    network_table,
    communities,
    community_column = "communities",
    weight_column = "weight") {
  # match TFs/TGs to their communities
  network_table$TG_module <- communities[, community_column][match(network_table$target, communities$gene)]
  network_table$TF_module <- communities[, community_column][match(network_table$regulator, communities$gene)]
  network_table$TG_module[is.na(network_table$TG_module)] <- network_table$TF_module[is.na(network_table$TG_module)]

  network_table$interaction <- "activation"
  network_table$interaction[network_table$corr < 0] <- "repression"

  flatDF <- network_table[, c("TF_module", "TG_module", weight_column, "interaction")]

  # make edge weight negative for repressive edges
  flatDF$edge_score <- flatDF[, weight_column]
  flatDF$edge_score[flatDF$interaction == "repression"] <- flatDF$edge_score[flatDF$interaction == "repression"] * -1

  # aggregate edges between modules by sum (repressive edges will just be negative)
  flatDF <- flatDF[, c("TF_module", "TG_module", "edge_score")]
  flatDF <- aggregate(. ~ TF_module + TG_module, data = flatDF, FUN = sum)

  # normalize edge weight by possible number of edges
  n_members <- data.frame(table(communities[, community_column]))
  flatDF$TG_members <- n_members$Freq[match(as.character(flatDF$TG_module), as.character(n_members$Var1))]
  flatDF$TF_members <- n_members$Freq[match(as.character(flatDF$TF_module), as.character(n_members$Var1))]

  flatDF$normalized_edge_score <- flatDF$edge_score / (flatDF$TF_members * flatDF$TG_members)
  flatDF$edge_length <- 1 - flatDF$normalized_edge_score

  flatDF[, c("TF_module", "TG_module", "edge_score", "normalized_edge_score", "edge_length")]
}


# network_table the result from running flatten_network
# module the module name at the end of the paths
find_paths_to <- function(network_table, module) {
  modnet_ig <- igraph::graph_from_data_frame(network_table, directed = TRUE)
  mods <- igraph::V(modnet_ig)$name
  mods <- mods[mods != module]

  res <- igraph::distances(modnet_ig, to = module, weights = igraph::E(modnet_ig)$edge_length, mode = "out")
  colnames(res) <- c("path_length")

  paths <- c()
  for (mod in rownames(res)) {
    path <- igraph::shortest_paths(modnet_ig, from = mod, to = module, mode = "out", weights = igraph::E(modnet_ig)$edge_length)$vpath[[1]]$name
    paths <- c(paths, paste(path, collapse = "--"))
  }

  res <- cbind(res, data.frame(path = paths))

  res$path[is.infinite(res$path_length)] <- NA

  res
}

#' @title rough_hierarchy
#'
#' @description
#'  returns rough roots in the network, rough roots selected as those connected to most number of nodes
#'
#' @param directed Whether the network is directed.
#' @inheritParams network_format
#'
#' @return list
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' rough_hierarchy(network_table)
rough_hierarchy <- function(
    network_table,
    abs_weight = TRUE,
    directed = TRUE) {
  # if (class(network_table) != "igraph") {
  if (!igraph::is.igraph(network_table)) {
    if (abs_weight) {
      network_table$weight <- abs(network_table$weight)
    }
    network_table <- igraph::graph_from_data_frame(
      network_table,
      directed = directed
    )
  }

  distances <- data.frame(
    igraph::distances(network_table, mode = "out")
  )
  distances$num_connected <- rowSums(distances != Inf)

  roots <- rownames(distances)[distances$num_connected == max(distances$num_connected)]

  return(
    list(
      roots = roots,
      num_paths = max(distances$num_connected)
    )
  )
}


# =================== Functions to perform shortest path analyses =====================

#' Function to return shortest path from 1 regulator to 1 target in a static network
#'
#' @param network_table a static network dataframe
#' @param regulator The starting gene.
#' @param target The end gene.
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#' @param compare_to_average if TRUE will compute normalized against average path length
#'
#' @return shortest path, distance, normalized distance, and action
#'
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' static_shortest_path(
#'   network_table,
#'   regulator = "g1",
#'   target = "g2"
#' )
static_shortest_path <- function(
    network_table,
    regulator,
    target,
    weight_column = "weight",
    compare_to_average = FALSE) {
  # compute relative edge lengths
  network_table$normalized_score <- normalization(
    network_table[, weight_column],
    method = "max"
  )
  network_table$edge_length <- 1 - network_table$normalized_score

  # find shortest paths
  ig <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "edge_length", "weight")],
    directed = TRUE
  )
  path <- igraph::shortest_paths(
    ig,
    from = regulator,
    to = target,
    mode = "out",
    output = "both",
    weights = igraph::E(ig)$edge_length
  )

  vpath <- path$vpath[[1]]$name
  epath <- path$epath[[1]]

  if (length(epath) == 0) {
    stop("No path")
  }

  # compute distance
  distance <- igraph::distances(
    ig,
    v = regulator,
    to = target,
    mode = "out",
    weights = igraph::E(ig)$edge_length
  )[1, 1]

  if (compare_to_average) {
    avg_path_length <- igraph::mean_distance(ig, directed = TRUE)
  }

  # dealing with repressive interactions
  # if path has even number of repressive interactions then regulator on = target turned on
  # else, if odd, regulator on = target off
  # nrep <- sum(path$epath[[1]]$corr < 0)

  if (compare_to_average) {
    list(
      path = vpath,
      distance = distance,
      distance_over_average = distance / avg_path_length,
      action = ifelse((sum(path$epath[[1]]$corr < 0) %% 2) == 0, 1, -1)
    )
  } else {
    list(
      path = vpath,
      distance = distance,
      action = ifelse((sum(path$epath[[1]]$corr < 0) %% 2) == 0, 1, -1)
    )
  }
}

# Dynamic version of ^^
# there is probably a better way of doing this that takes in to account the dynamic aspect of the network
# for now, this function will merge networks and run shortest path on merged network..

#' Function to return shortest path from 1 regulator to 1 target in a dynamic network
#'
#' @param network_table a dyanmic network
#' @param regulator the starting regulator
#' @param target the end regulator
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#' @param compare_to_average if TRUE will compute normalized against average path length
#'
#' @return shortest path, distance, normalized distance, and action
#'
#' @export
dynamic_shortest_path <- function(
    network_table,
    regulator,
    target,
    weight_column = "weight",
    compare_to_average = FALSE) {
  network_table <- do.call("rbind", network_table)
  network_table <- network_table[!duplicated(network_table[, c("regulator", "target")]), ]

  return(
    static_shortest_path(
      network_table,
      regulator,
      target,
      weight_column,
      compare_to_average
    )
  )
}

# quick function to loop through ^^ for multiple TFs and targets. Returns a data frame
#
# network_table dynamic network (list)
# from a vector of TFs
# to a vector of targets

#' Function to return shortest path from multiple TFs to multiple targets in a dynamic network
#'
#' @param network_table a dyanmic network
#' @param from the starting TFs
#' @param targets the end TFs
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#'
#' @return dataframe with shortest path, distance, normalized distance, and action
#'
#' @export
dynamic_shortest_path_multiple <- function(
    network_table,
    regulators,
    targets,
    weight_column = "weight") {
  res <- data.frame(
    regulator = character(),
    target = character(),
    path = character(),
    distance = numeric(),
    action = numeric()
  )

  for (regulator in regulators) {
    for (target in targets) {
      tryCatch(
        {
          path <- dynamic_shortest_path(
            network_table,
            regulator = regulator,
            target = target,
            weight_column = weight_column,
            compare_to_average = FALSE
          )
          res <- rbind(
            res,
            data.frame(
              regulator = regulator,
              target = target,
              path = paste(path[["path"]], collapse = "--"),
              distance = path[["distance"]],
              action = path[["action"]]
            )
          )
        },
        error = function(e) {}
      )
    }
  }

  network_table <- do.call("rbind", network_table)
  network_table$normalized_score <- normalization(
    network_table[, weight_column],
    method = "max"
  )
  network_table$edge_length <- 1 - network_table$normalized_score
  ig <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "edge_length", "corr")],
    directed = TRUE
  )

  avg_path_length <- igraph::mean_distance(ig, directed = TRUE)
  res$distance_over_average <- res$distance / avg_path_length

  return(res)
}

#' Adds an extra column to the result of dynamic_shortest_path_multiple that predicts overall action based on correlation between "from" and "to"
#'
#' @param spDF the result of running dynamic_shortest_path_multiple
#' @param matrix corresponding genes-by-cells expression matrix
#'
#' @return shortest paths dataframe with added action_by_corr column
#'
#' @export
cor_and_add_action <- function(
    spDF,
    matrix) {
  genes <- union(spDF$from, spDF$to)
  matrix <- matrix[genes, ]

  corrs <- apply(spDF[, c("from", "to")], 1, function(x) {
    cor(matrix[x[1], ], matrix[x[2], ])
  })
  spDF$action_by_corr <- corrs
  spDF$action_by_corr <- ifelse(spDF$action_by_corr < 0, -1, 1)

  spDF
}

# =================== Functions to perform reachability analyses =====================

# a very basic reachability function....
# max dist numeric value, otherwise "less_than_mean"
# from is a either a character or vector of characters
static_reachability <- function(
    network_table,
    from,
    max_dist = "less_than_mean",
    weight_column = "weight",
    tfs = NULL,
    tf_only = FALSE) {
  # compute relative edge lengths
  network_table$normalized_score <- normalization(
    network_table[, weight_column],
    method = "max"
  )
  network_table$edge_length <- 1 - network_table$normalized_score

  if (tf_only) {
    if (is.null(tfs)) {
      stop("Supply TFs.")
    } else {
      network_table <- network_table[network_table$target %in% tfs, ]
    }
  }

  ig <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "edge_length", "corr")],
    directed = TRUE
  )

  if (max_dist == "less_than_mean") {
    avg_path_length <- igraph::mean_distance(ig, directed = TRUE)
    max_dist <- avg_path_length
  }

  from <- intersect(from, union(network_table$regulator, network_table$target))

  # compute distances
  distance_mat <- igraph::distances(
    ig,
    v = from,
    weights = igraph::E(ig)$edge_length,
    mode = "out"
  )

  reachable <- lapply(1:nrow(distance_mat), function(x) {
    t <- colnames(distance_mat)[(distance_mat[x, ] < max_dist)]
    t <- t[t != rownames(distance_mat)[x]]
    t
  })
  names(reachable) <- rownames(distance_mat)

  reachable
}

# the dynamic version of ^^ .. again, initial simplified version
dynamic_reachability <- function(
    network_table,
    from,
    max_dist = "less_than_mean",
    weight_column = "weight",
    tfs = NULL,
    tf_only = FALSE) {
  # merge
  network_table <- do.call("rbind", network_table)
  network_table <- network_table[!duplicated(network_table[, c("regulator", "target")]), ]

  res <- static_reachability(
    network_table,
    from = from,
    max_dist = max_dist,
    weight_column = weight_column,
    tfs = tfs,
    tf_only = tf_only
  )
  res
}


# =================== Useful plotting functions =====================
# function to plot heatmap with pre-split matrix
# expList list of expression matrices

#' Useful plotting function to plot heatmap of module expression across time with pre-split matrix
#'
#' @param expList list of expression matrices
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param toScale whether or not to scale expression
#' @param limits limits on expression
#' @param smooth whether or not to smooth expression across pseudotime for cleaner plotting
#' @param order_by name in expList that is used to order rows in the heatmap
#' @param thresh_on threshold expression is considered on, used in ordering the rows
#' @param fontSize heatmap font size
#' @param anno_colors annotation colors
#'
#' @return pheatmap
#'
#' @export
heatmap_by_treatment_group <- function(
    expList,
    meta_data,
    pseudotime_column = "pseudotime",
    toScale = T,
    limits = c(0, 5),
    smooth = TRUE,
    order_by = "WAG",
    thresh_on = 0.02,
    fontSize = 8,
    anno_colors = NA) {
  # if (class(expList) != "list") {
  if (!is.list(expList)) {
    expList <- list(A = expList)
  }
  expList <- lapply(expList, function(x) {
    st <- meta_data[colnames(x), ]
    st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
    x[, rownames(st)]
  })
  if (smooth) {
    expList <- lapply(expList, function(x) {
      st <- data.frame(
        cell_name = rownames(meta_data[colnames(x), ]),
        pseudotime = meta_data[colnames(x), pseudotime_column]
      )
      rownames(st) <- st$cell_name
      expression_ksmooth(x, st, bandwith = 0.05)
    })
  }

  expList <- lapply(expList, function(x) {
    x[rowSums(is.na(x)) == 0, ]
  })

  # For this function, order rows by time of hitting some expression threshold
  peak_time <- apply(expList[[order_by]], 1, function(x) {
    ifelse(any(x > thresh_on), mean(which(x > thresh_on)[1:10]), length(x) + 1)
  })
  peak_time[is.na(peak_time)] <- ncol(expList[[order_by]]) + 1

  # peak_time = apply(expList[[order_by]], 1, which.max)
  genes_ordered <- names(sort(peak_time))

  matrix <- do.call(cbind, expList)
  meta_data <- meta_data[colnames(matrix), ]

  meta_data$cell_name <- rownames(meta_data)
  if ("epoch" %in% colnames(meta_data)) {
    col_ann <- meta_data[, c("treatment", pseudotime_column, "epoch")]
  } else {
    col_ann <- meta_data[, c("treatment", pseudotime_column)]
  }

  value <- matrix[genes_ordered, ]
  if (toScale) {
    if (class(value)[1] != "matrix") {
      value <- t(scale(Matrix::t(value)))
    } else {
      value <- t(scale(t(value)))
    }
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

  # val_col <- grDevices::colorRampPalette(
  #   rev(RColorBrewer::brewer.pal(n = 12, name = "Spectral"))
  # )(25)

  # find gaps rows and columns
  gaps <- sapply(expList, ncol)
  gaps <- cumsum(gaps)
  gaps <- gaps[-length(gaps)]

  pheatmap::pheatmap(
    value,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = F,
    annotation_colors = anno_colors,
    annotation_col = col_ann,
    border_color = NA,
    gaps_col = gaps,
    annotation_names_row = F,
    fontsize = fontSize
  )
}

#' Useful plotting function to plot heatmap with pre-split matrix
#'
#' @param expList list of expression matrices
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param toScale whether or not to scale expression
#' @param limits limits on expression
#' @param smooth whether or not to smooth expression across pseudotime for cleaner plotting
#' @param fontSize heatmap font size
#' @param anno_colors annotation colors
#'
#' @return pheatmap
#'
#' @export
plot_heatmap_by_treatment <- function(
    expList,
    meta_data,
    pseudotime_column = "latent_time",
    toScale = T,
    limits = c(0, 5),
    smooth = TRUE,
    anno_colors = NA,
    fontSize = 8) {
  # if (class(expList) != "list") {
  if (!is.list(expList)) {
    expList <- list(A = expList)
  }

  expList <- lapply(expList, function(x) {
    st <- meta_data[colnames(x), ]
    st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
    x[, rownames(st)]
  })
  if (smooth) {
    expList <- lapply(expList, function(x) {
      st <- data.frame(cell_name = rownames(meta_data[colnames(x), ]), pseudotime = meta_data[colnames(x), pseudotime_column])
      rownames(st) <- st$cell_name
      expression_ksmooth(x, st, bandwith = 0.03)
    })
  }

  # For this function, order rows by time of peak expression in WAG (first element in list)
  # chunked by community epochs
  rows <- data.frame(subnet = rownames(expList[[1]]), network = unlist(lapply(strsplit(rownames(expList[[1]]), "_"), "[[", 1)))
  genes_ordered <- c()
  for (n in unique(rows$network)) {
    peak_time <- apply(expList[[1]][startsWith(rownames(expList[[1]]), n), ], 1, function(x) {
      mean(order(x, decreasing = TRUE)[1:30])
    })
    genes_ordered <- c(genes_ordered, names(sort(peak_time)))
  }

  matrix <- do.call(cbind, expList)
  meta_data <- meta_data[colnames(matrix), ]

  meta_data$cell_name <- rownames(meta_data)
  col_ann <- meta_data[, c("treatment", pseudotime_column)]

  row_ann <- data.frame(subnet = rownames(matrix), network = unlist(lapply(strsplit(rownames(matrix), "_"), "[[", 1)))
  rownames(row_ann) <- row_ann$subnet
  row_ann$subnet <- NULL
  row_ann$network <- factor(row_ann$network, levels = unique(row_ann$network))

  value <- matrix[genes_ordered, ]
  if (toScale) {
    if (class(value)[1] != "matrix") {
      value <- t(scale(Matrix::t(value)))
    } else {
      value <- t(scale(t(value)))
    }
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

  # val_col <- grDevices::colorRampPalette(
  #   rev(RColorBrewer::brewer.pal(n = 12, name = "Spectral"))
  # )(25)

  # find gaps rows and columns
  gaps <- sapply(expList, ncol)
  gaps <- cumsum(gaps)
  gaps <- gaps[-length(gaps)]

  rowgaps <- table(row_ann$network)
  rowgaps <- cumsum(rowgaps)
  rowgaps <- rowgaps[-length(rowgaps)]

  pheatmap::pheatmap(
    value,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = F,
    annotation_colors = anno_colors,
    annotation_col = col_ann,
    annotation_row = row_ann,
    border_color = NA,
    gaps_col = gaps,
    gaps_row = rowgaps,
    annotation_names_row = F,
    fontsize = fontSize
  )
}

plot_reachability_community_coverage <- function(
    reachableList,
    communities) {
  # get rid of duplicates in reachableList, in case there are any
  reachableList <- lapply(reachableList, function(x) {
    x[!duplicated(x)]
  })

  res <- data.frame(matrix(ncol = length(reachableList), nrow = length(communities)))
  colnames(res) <- names(reachableList)
  rownames(res) <- names(communities)

  for (m in names(communities)) {
    pct_covered <- sapply(reachableList, function(x) {
      sum(x %in% communities[[m]]) / length(communities[[m]])
    })
    res[m, ] <- pct_covered
  }

  res
}

#' Function that orders genes based on peak expression
#'
#' @param matrix expression matrix
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param smooth whether or not to smooth expression across pseudotime
#'
#' @return ordered genes
#'
#' @export
order_genes <- function(
    matrix,
    meta_data,
    pseudotime_column,
    smooth = TRUE) {
  st <- meta_data[colnames(matrix), ]
  st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
  matrix <- matrix[, rownames(st)]

  if (smooth) {
    st <- data.frame(
      cell_name = rownames(st),
      pseudotime = st[, pseudotime_column]
    )
    rownames(st) <- st$cell_name
    matrix <- expression_ksmooth(matrix, st, bandwith = 0.05)
  }

  peak_time <- apply(matrix, 1, which.max)
  genes_ordered <- names(sort(peak_time))

  return(genes_ordered)
}
