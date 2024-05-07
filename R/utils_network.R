#' Function to compute page rank of TF+target networks
#'
#' @param networks_list result of dynamic network reconstruction
#' @param regulators list of regulators
#' @param targets list of targets
#' @param directed if network is directed or not
#'
#' @return list of each network, ranked by `page_rank`
#' @export
compute_pagerank <- function(
    networks_list,
    regulators = NULL,
    targets = NULL,
    directed = FALSE) {
  ranked_networks <- sapply(
    names(networks_list), function(x) NULL
  )
  for (net in names(networks_list)) {
    ranked_networks[[net]] <- calculate.gene.rank(
      networks_list[[net]],
      regulators = regulators,
      targets = targets
    )
  }

  return(ranked_networks)
}

#' Function to compute betweenness and degree
#'
#' @inheritParams compute_pagerank
#'
#' @return list of each network, ranked by `betweenness*degree`
#' @export
compute_betweenness_degree <- function(
    networks_list,
    regulators = NULL,
    targets = NULL,
    directed = FALSE) {
  ranked_networks <- sapply(names(networks_list), function(x) NULL)
  for (net_name in names(networks_list)) {
    df <- networks_list[[net_name]]
    df <- net.format(
      df,
      regulators = regulators,
      targets = targets
    )

    net <- igraph::graph_from_data_frame(df, directed = directed)
    b <- igraph::betweenness(net, directed = directed, normalized = TRUE)
    b <- as.data.frame(b)

    degree <- igraph::degree(net, mode = "all", normalized = TRUE)
    degree <- data.frame(degree[rownames(b)])

    b <- cbind(b, degree)
    colnames(b) <- c("betweenness", "degree")
    b$temp <- b$betweenness * b$degree

    b <- b[order(b$temp, decreasing = TRUE), ]
    b$temp <- NULL

    b$gene <- rownames(b)
    b <- b[, c("gene", "betweenness", "degree")]

    b$is_regulator <- FALSE
    b$is_regulator[b$gene %in% unique(df$regulator)] <- TRUE

    ranked_networks[[net_name]] <- b
  }

  return(ranked_networks)
}

#' @title add_interaction_type
#' @description
#'  Add interaction type to static or dynamic network
#'
#' @param network a network table/list
#'
#' @return a network table/list with network_add interaction type
#'
#' @export
add_interaction_type <- function(network) {
  # if (class(network) == "data.frame") {
  if (is.data.frame(network)) {
    if (!("weight" %in% names(network))) {
      stop(
        "Missing 'weight' information.
        Run `inferCSN` first."
      )
    } else {
      network_add <- network
      network_add$interaction <- ifelse(
        network_add$weight < 0, "Repression", "Activation"
      )
    }
  # } else if (class(network) == "list") {
  } else if (is.list(network)) {
    keep <- sapply(network, function(x) {
      ("weight" %in% names(x))
    })
    if (any(keep == FALSE)) {
      stop(
        "At least one network is missing 'weight' information.
        Run `inferCSN` first."
      )
    } else {
      network_add <- lapply(
        network,
        transform,
        interaction = ifelse(weight < 0, "Repression", "Activation")
      )
    }
  } else {
    stop("Invalid network.")
  }

  return(network_add)
}

#' @title expression_ksmooth
#'
#' @description
#'  smooths expression across cells in path.
#'
#' @param matrix expression matrix.
#' @param meta_data meta_data includes "cell_name", "pseudotime", "group".
#' @param pseudotime_column pseudotime_column.
#' @param bandwith bandwith for kernel smoother.
#'
#' @return smoothed matrix
#'
#' @export
expression_ksmooth <- function(
    matrix,
    meta_data,
    pseudotime_column = "pseudotime",
    bandwith = 0.25) {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  cells_used <- intersect(meta_data$cells, colnames(matrix))
  meta_data <- meta_data[cells_used, ]
  bandwith <- min(
    bandwith,
    (max(meta_data$pseudotime) - min(meta_data$pseudotime)) / 10
  )

  pseudotime <- meta_data$pseudotime
  matrix <- matrix[, cells_used]

  matrix_new <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  for (i in seq_len(nrow(matrix))) {
    matrix_new[i, ] <- stats::ksmooth(
      pseudotime,
      matrix[i, ],
      kernel = "normal",
      bandwidth = bandwith,
      x.points = pseudotime
    )$y
  }

  rownames(matrix_new) <- rownames(matrix)
  colnames(matrix_new) <- colnames(matrix)

  return(matrix_new)
}

#' @title pre_pseudotime_matrix
#'
#' @description
#'  Calculate and return a smoothed pseudotime matrix for the given gene list
#'
#' @param matrix matrix
#' @param gene_list A vector of gene names
#'
#' @return A smoothed pseudotime matrix for the given gene list
#'
#' @export
pre_pseudotime_matrix <- function(
    matrix,
    gene_list = NULL) {
  if (is.null(gene_list)) {
    gene_list <- rownames(matrix)
  }
  # choose assays type
  pt_matrix <- matrix

  pt_matrix <- t(apply(
    pt_matrix, 1, function(x) {
      stats::smooth.spline(x, df = 3)$y
    }
  ))
  pt_matrix <- t(apply(
    pt_matrix, 1, function(x) {
      (x - mean(x)) / stats::sd(x)
    }
  ))
  rownames(pt_matrix) <- gene_list
  colnames(pt_matrix) <- colnames(matrix)

  return(pt_matrix)
}

#' @title find_communities
#'
#' @description
#'  Find communities in a static or dynamic network
#'
#' @param weight_table weight_table
#' @param use_weights whether or not to use edge weights (for weighted graphs)
#' @param directed Default set to `FALSE`.
#'  If directed = TRUE, remember to flip target and regulator in diffnet dfs
#'
#' @return community assignments of nodes in the dynamic network
#'
#' @export
find_communities <- function(
    weight_table,
    use_weights = FALSE,
    directed = FALSE) {
  # if (class(weight_table) == "data.frame") {
  if (is.data.frame(weight_table)) {
    graphs <- igraph::graph_from_data_frame(weight_table, directed = directed)

    if (use_weights) {
      weights <- igraph::edge_attr(graphs, "weight")
    } else {
      weights <- NA
    }

    communities <- as.data.frame(
      as.table(
        igraph::membership(
          igraph::cluster_louvain(graphs, weights = weights)
        )
      )
    )
    colnames(communities) <- c("gene", "communities")
  # } else if (class(weight_table) == "list") {
  } else if (is.list(weight_table)) {
    if (use_weights) {
      graphs <- lapply(
        weight_table, function(x) {
          igraph::graph_from_data_frame(
            x[, c("regulator", "target", "weight")],
            directed = directed
          )
        }
      )
      communities <- lapply(
        graphs, function(x) {
          weights <- igraph::edge_attr(x, "weight")
          c <- as.data.frame(
            as.table(
              igraph::membership(
                igraph::cluster_louvain(x, weights = weights)
              )
            )
          )
          colnames(c) <- c("gene", "communities")
          c
        }
      )
    } else {
      graphs <- lapply(
        weight_table, function(x) {
          igraph::graph_from_data_frame(x, directed = FALSE)
        }
      )
      communities <- lapply(graphs, function(x) {
        c <- as.data.frame(
          as.table(
            igraph::membership(
              igraph::cluster_louvain(x, weights = NA)
            )
          )
        )
        colnames(c) <- c("gene", "communities")
        c
      })
    }
  } else {
    stop("Invalid network.")
  }

  return(communities)
}


#' @title compileDynGenes
#'
#' @description
#'  compiles results from findDynGenes run on separate paths into a single result
#'
#' @param dynReslist list of results from findDynGenes
#'
#' @return a new dynRes a.k.a. the compiled result of findDynGenes
#'
#' @export
compileDynGenes <- function(dynReslist) {
  gene_df <- as.data.frame(
    sapply(dynReslist, function(x) x[["genes"]])
  )
  gene_df$min_pval <- apply(gene_df, 1, min)
  genes <- c(gene_df$min_pval)
  names(genes) <- rownames(gene_df)

  cell_list <- lapply(dynReslist, function(x) x[["cells"]])
  cell_df <- as.data.frame(do.call(rbind, cell_list))
  cell_df <- cell_df[!duplicated(cell_df), ]
  cell_df <- cell_df[order(cell_df$pseudotime, decreasing = FALSE), ]

  list(genes = genes, cells = cell_df)
}

#' @title compile_epochs
#'
#' @description
#'  compile networks based on matching network label in result of running assign_network.
#'
#' @param networks_list list of results of running assign_network on each path.
#'
#' @return a compiled list of network genes
#' @export
compile_epochs <- function(networks_list) {
  compiled_networks <- list()

  for (path in 1:length(networks_list)) {
    networks <- networks_list[[path]]
    networks$mean_expression <- NULL
    for (network in names(networks)) {
      if (network %in% names(compiled_networks)) {
        compiled_networks[[network]] <- union(
          networks[[network]],
          compiled_networks[[network]]
        )
      } else {
        compiled_networks[[network]] <- networks[[network]]
      }
    }
  }

  means <- data.frame(
    gene = character(),
    network = character(),
    mean_expression = numeric()
  )
  for (path in 1:length(networks_list)) {
    mean_exp <- networks_list[[path]]$mean_expression
    means <- rbind(means, mean_exp)
  }
  compiled_means <- aggregate(
    means$mean_expression,
    by = list(gene = means$gene, network = means$network),
    data = means,
    FUN = mean
  )
  colnames(compiled_means) <- c("gene", "network", "mean_expression")

  compiled_networks[["mean_expression"]] <- compiled_means

  return(compiled_networks)
}
