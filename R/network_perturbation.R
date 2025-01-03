#' Set Network object
#' @param object CSNObject
#' @param network Network object to set
#' @param name Character, name of the network (optional)
#' @return Updated CSNObject
#' @export
SetNetwork <- function(object, network, name = NULL) {
  if (is.null(name)) {
    name <- object@active_network
  }

  if (is.null(name)) {
    stop("No active network specified")
  }

  object@networks[[name]] <- network
  object@active_network <- name
  return(object)
}

#' Predict gene expression under TF perturbation
#' @param object A Network or CSNObject object
#' @param perturb_tfs Named vector of TF perturbation states
#' @param use_weight Logical, whether to use network weights instead of coefficients
#' @param n_iter Number of iterations for perturbation simulation
#' @param scale_factor Scale factor for perturbation effect (0-1)
#' @param verbose Whether to show progress messages
#' @return Updated object with perturbation predictions
#' @export
setGeneric("predictPerturbation", function(object,
                                           perturb_tfs,
                                           use_weight = FALSE,
                                           n_iter = 5,
                                           scale_factor = 0.1,
                                           verbose = TRUE) {
  standardGeneric("predictPerturbation")
})

#' @rdname predictPerturbation
#' @export
setMethod(
  "predictPerturbation",
  signature(object = "Network"),
  function(object, perturb_tfs, use_weight = FALSE, n_iter = 3,
           scale_factor = 0.1, verbose = TRUE) {
    if (is.null(names(perturb_tfs))) {
      stop("perturb_tfs must be a named vector")
    }
    if (!all(names(perturb_tfs) %in% object@regulators)) {
      stop("Some perturbed TFs are not in the network regulators.")
    }
    if (scale_factor <= 0 || scale_factor > 1) {
      stop("'scale_factor' must be between 0 and 1.")
    }
    if (n_iter <= 0 || !is.integer(as.integer(n_iter))) {
      stop("'n_iter' must be a positive integer.")
    }

    tf_names <- names(perturb_tfs)
    perturb_values <- as.numeric(unname(perturb_tfs))
    names(perturb_values) <- tf_names

    if (any(is.na(perturb_values))) {
      stop("All perturbation values must be numeric")
    }

    orig_expr <- as.matrix(object@data)
    pert_expr <- orig_expr

    if (use_weight) {
      net_df <- object@network
      B <- Matrix::sparseMatrix(
        i = match(net_df$target, colnames(orig_expr)),
        j = match(net_df$regulator, colnames(orig_expr)),
        x = net_df$weight,
        dims = c(ncol(orig_expr), ncol(orig_expr))
      )
    } else {
      coef_df <- object@coefficients
      B <- Matrix::sparseMatrix(
        i = match(coef_df$target, colnames(orig_expr)),
        j = match(coef_df$tf, colnames(orig_expr)),
        x = coef_df$coefficient,
        dims = c(ncol(orig_expr), ncol(orig_expr))
      )
    }

    delta_input <- matrix(0, nrow = nrow(orig_expr), ncol = ncol(orig_expr))
    colnames(delta_input) <- colnames(orig_expr)

    for (i in seq_along(tf_names)) {
      tf <- tf_names[i]
      target_value <- perturb_values[i]
      delta_input[, tf] <- target_value - orig_expr[, tf]
    }

    current_delta <- delta_input

    for (i in 2:n_iter) {
      if (verbose) message(sprintf("Iteration %d/%d", i, n_iter))

      current_delta <- as.matrix(current_delta %*% B)

      current_expr <- orig_expr + scale_factor * current_delta

      current_expr[current_expr < 0] <- 0
    }

    tf_states <- data.frame(
      tf = tf_names,
      state = perturb_values,
      stringsAsFactors = FALSE,
      row.names = NULL,
      check.names = FALSE
    )

    pert_obj <- new("Perturbation",
      original = as.matrix(orig_expr),
      perturbed = as.matrix(current_expr),
      tf_states = tf_states,
      params = list(
        method = if (use_weight) "network_weight" else "coefficient",
        n_iterations = n_iter,
        scale_factor = scale_factor,
        final_delta = current_delta
      )
    )

    object@perturbation <- pert_obj
    return(object)
  }
)

#' @rdname predictPerturbation
#' @export
setMethod(
  "predictPerturbation", signature(object = "CSNObject"),
  function(object, perturb_tfs, use_weight = FALSE) {
    net <- GetNetwork(object)
    net <- predictPerturbation(net, perturb_tfs, use_weight)
    object@networks[[object@active_network]] <- net
    return(object)
  }
)

#' Embed perturbation results
#' @param object A Network or CSNObject object
#' @param reduction_method Character, dimensionality reduction method ("umap", "tsne", "pca")
#' @param dims Number of dimensions to compute
#' @param scale Whether to scale the data before reduction
#' @param cores Number of cores for parallel processing
#' @param seed Random seed for reproducibility
#' @return Updated object with embedded perturbation results
#' @export
setGeneric("embedPerturbation", function(object,
                                         reduction_method = "umap",
                                         dims = 2,
                                         scale = TRUE,
                                         cores = 1,
                                         seed = 1) {
  standardGeneric("embedPerturbation")
})

#' @rdname embedPerturbation
#' @export
setMethod(
  "embedPerturbation",
  signature(object = "Network"),
  function(object, reduction_method = "umap", dims = 2, scale = TRUE,
           cores = 1, seed = 1) {
    pert <- object@perturbation

    orig_data <- pert@original
    pert_data <- pert@perturbed
    delta_expr <- pert_data - orig_data
    combined_expr <- rbind(orig_data, pert_data)

    if (is.null(rownames(combined_expr))) {
      rownames(combined_expr) <- c(
        paste0("orig_", seq_len(nrow(orig_data))),
        paste0("pert_", seq_len(nrow(pert_data)))
      )
    }

    combined_expr <- preprocess_data(combined_expr, scale)

    embedding <- tryCatch(
      {
        do_reduction(combined_expr, reduction_method, dims, cores, seed)
      },
      error = function(e) {
        warning(sprintf("%s failed, falling back to PCA", reduction_method))
        do_reduction(combined_expr, "pca", dims, cores, seed)
      }
    )

    n_orig <- nrow(orig_data)
    orig_embedding <- embedding[1:n_orig, , drop = FALSE]
    pert_embedding <- embedding[(n_orig + 1):nrow(embedding), , drop = FALSE]

    trajectory_vectors <- calculate_trajectory(
      embedding = orig_embedding,
      delta_expr = delta_expr,
      n_neighbors = min(15, nrow(orig_embedding) - 1)
    )

    pert_trajectories <- matrix(0, nrow = nrow(pert_embedding), ncol = ncol(embedding))
    full_trajectories <- rbind(trajectory_vectors, pert_trajectories)

    pert@embedding <- embedding
    pert@params$reduction <- list(
      method = reduction_method,
      dims = dims,
      scale = scale,
      seed = seed,
      data_shape = dim(combined_expr),
      var_features = colnames(combined_expr),
      n_orig = n_orig,
      n_pert = nrow(pert_data),
      trajectory_vectors = full_trajectories
    )

    object@perturbation <- pert
    return(object)
  }
)

preprocess_data <- function(data, scale) {
  if (any(is.na(data))) data[is.na(data)] <- 0

  var_cols <- apply(data, 2, var) > 1e-10
  if (!any(var_cols)) stop("No variable columns found in the data")
  data <- data[, var_cols]

  if (scale) {
    col_means <- colMeans(data)
    col_sds <- apply(data, 2, sd)
    col_sds[col_sds < 1e-10] <- 1
    data <- scale(data, center = col_means, scale = col_sds)
  }

  return(data)
}

do_reduction <- function(data, method, dims, cores, seed) {
  if (dims > ncol(data)) {
    stop("Number of dimensions exceeds data feature count.")
  }
  if (dims < 1) {
    stop("Number of dimensions must be positive.")
  }

  set.seed(seed)
  result <- switch(method,
    "umap" = uwot::umap(
      data,
      n_components = dims,
      n_threads = cores,
      seed = seed,
      min_dist = 0.1,
      n_neighbors = min(15, nrow(data) - 1)
    ),
    "tsne" = Rtsne::Rtsne(
      data,
      dims = dims,
      num_threads = cores,
      perplexity = min(30, nrow(data) / 4),
      check_duplicates = FALSE,
      pca = TRUE
    )$Y,
    "pca" = {
      pca_result <- stats::prcomp(data, scale. = TRUE)
      pca_result$x[, seq_len(dims), drop = FALSE]
    },
    stop("Unsupported reduction method: ", method)
  )

  colnames(result) <- paste0("Dim", seq_len(dims))
  rownames(result) <- rownames(data)

  return(result)
}

#' Visualize perturbation results
#' @param object A Network or CSNObject object
#' @param label_points Logical, whether to add point labels
#' @param point_size Numeric, size of points
#' @param max_overlaps Maximum number of overlapping labels to allow
#' @param alpha Point transparency
#' @param show_trajectory Logical, whether to show trajectory arrows
#' @param arrow_size Size of trajectory arrows
#' @param grid_n Number of grid points for trajectory arrows
#' @return ggplot object
#' @export
setGeneric("plotPerturbation", function(object,
                                        label_points = FALSE,
                                        point_size = 1,
                                        max_overlaps = 10,
                                        alpha = 0.6,
                                        show_trajectory = TRUE,
                                        arrow_size = 0.5,
                                        grid_n = 25) {
  standardGeneric("plotPerturbation")
})

#' @rdname plotPerturbation
#' @export
setMethod(
  "plotPerturbation",
  signature(object = "Network"),
  function(object, label_points = FALSE, point_size = 1,
           max_overlaps = 10, alpha = 0.6,
           show_trajectory = TRUE, arrow_size = 0.5,
           grid_n = 25) {
    pert <- object@perturbation

    if (is.null(pert@embedding)) {
      message("No embedding found. Running embedPerturbation...")
      object <- embedPerturbation(object)
      pert <- object@perturbation
    }

    n_orig <- nrow(pert@original)
    n_pert <- nrow(pert@perturbed)

    plot_df <- data.frame(
      DR1 = pert@embedding[, 1],
      DR2 = pert@embedding[, 2],
      Type = factor(
        c(rep("Original", n_orig), rep("Perturbed", n_pert)),
        levels = c("Original", "Perturbed")
      )
    )

    p <- ggplot(plot_df, aes(x = DR1, y = DR2, color = Type)) +
      geom_point(alpha = alpha, size = point_size) +
      scale_color_manual(values = c(
        "Original" = "#66CC33",
        "Perturbed" = "#9933CC"
      )) +
      theme_minimal()

    if (show_trajectory && !is.null(pert@params$reduction$trajectory_vectors)) {
      traj_vectors <- pert@params$reduction$trajectory_vectors
      arrow_df <- data.frame(
        x = plot_df$DR1,
        y = plot_df$DR2,
        xend = plot_df$DR1 + traj_vectors[, 1],
        yend = plot_df$DR2 + traj_vectors[, 2]
      )

      p <- p + geom_segment(
        data = arrow_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(arrow_size, "cm")),
        color = "grey30",
        alpha = 0.5
      )
    }

    p <- p +
      labs(
        title = "Perturbation Analysis",
        subtitle = sprintf(
          "Method: %s (%d iterations)",
          pert@params$method,
          pert@params$n_iterations
        ),
        x = "Dimension 1",
        y = "Dimension 2"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12, color = "grey40"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA)
      )

    return(p)
  }
)

#' @rdname plotPerturbation
#' @export
setMethod(
  "plotPerturbation", signature(object = "CSNObject"),
  function(object, label_points = FALSE, point_size = 1,
           max_overlaps = 10, alpha = 0.6) {
    net <- GetNetwork(object)
    plotPerturbation(net, label_points, point_size, max_overlaps, alpha)
  }
)

#' Plot perturbation trajectory
#' @param object Network or CSNObject
#' @param genes Genes to plot
#' @param show_heatmap Logical, whether to show expression heatmap
#' @param n_top_genes Number of top variable genes to show in heatmap
#' @param heatmap_style Style of heatmap visualization ("combined" or "separate")
#' @param ... Additional arguments
#' @return ggplot object or list of plots
#' @export
setGeneric("plotPerturbationTrajectory", function(object,
                                                  genes = NULL,
                                                  show_heatmap = FALSE,
                                                  n_top_genes = 50,
                                                  heatmap_style = "combined",
                                                  ...) {
  standardGeneric("plotPerturbationTrajectory")
})

#' @rdname plotPerturbationTrajectory
#' @export
setMethod(
  "plotPerturbationTrajectory",
  signature(object = "Network"),
  function(object, genes = NULL, show_heatmap = FALSE, n_top_genes = 50, 
           heatmap_style = "combined", ...) {
    pert <- object@perturbation
    
    if (is.null(genes)) {
      genes <- pert@tf_states$tf
    }
    
    orig_expr <- pert@original
    pert_expr <- pert@perturbed
    
    delta_expr <- pert_expr - orig_expr
    
    cv <- apply(delta_expr, 2, function(x) sd(x) / mean(abs(x)))
    
    top_genes <- names(sort(cv, decreasing = TRUE))[1:n_top_genes]
    
    if (heatmap_style == "separate") {
      heatmap_data <- cbind(
        orig_expr[, top_genes],
        pert_expr[, top_genes]
      )
      
      col_names <- c(
        paste0("orig_", top_genes),
        paste0("pert_", top_genes)
      )
      colnames(heatmap_data) <- col_names
      
      column_anno <- data.frame(
        State = c(
          rep("Original", length(top_genes)),
          rep("Perturbed", length(top_genes))
        ),
        Gene = c(
          top_genes,
          top_genes
        ),
        row.names = col_names
      )
      
      h <- pheatmap::pheatmap(
        heatmap_data,
        scale = "row",
        clustering_method = "ward.D2",
        show_rownames = FALSE,
        show_colnames = TRUE,
        annotation_col = column_anno,
        gaps_col = length(top_genes),
        main = "Expression Changes Heatmap (Original vs Perturbed)",
        annotation_colors = list(
          State = c(
            Original = "#66CC33",
            Perturbed = "#9933CC"
          )
        )
      )
    } else {
      heatmap_data <- cbind(
        orig_expr[, top_genes],
        pert_expr[, top_genes]
      )
      
      column_anno <- data.frame(
        State = rep(c("Original", "Perturbed"), each = length(top_genes)),
        Gene = rep(top_genes, 2)
      )
      
      col_names <- paste(
        rep(c("orig", "pert"), each = length(top_genes)),
        rep(top_genes, 2),
        sep = "_"
      )
      colnames(heatmap_data) <- col_names
      rownames(column_anno) <- col_names
      
      h <- pheatmap::pheatmap(
        heatmap_data,
        scale = "row",
        clustering_method = "ward.D2",
        show_rownames = FALSE,
        show_colnames = FALSE,
        annotation_col = column_anno,
        main = "Expression Changes Heatmap",
        annotation_colors = list(
          State = c(
            Original = "#66CC33",
            Perturbed = "#9933CC"
          )
        )
      )
    }
    
    traj_data <- data.frame(
      Iteration = rep(c(1, 2), each = nrow(orig_expr) * length(genes)),
      Gene = rep(rep(genes, each = nrow(orig_expr)), 2),
      Expression = c(
        as.vector(orig_expr[, genes]),
        as.vector(pert_expr[, genes])
      ),
      Type = rep(c("Original", "Perturbed"), each = nrow(orig_expr) * length(genes))
    )
    
    p <- ggplot(traj_data, aes(x = Iteration, y = Expression, color = Type)) +
      geom_boxplot(alpha = 0.6) +
      facet_wrap(~Gene, scales = "free_y") +
      scale_color_manual(values = c("Original" = "#66CC33", "Perturbed" = "#9933CC")) +
      theme_minimal() +
      labs(
        title = "Gene Expression Trajectories",
        x = "State",
        y = "Expression Level"
      )
    
    if (show_heatmap) {
      return(list(trajectory = p, heatmap = h))
    }
    
    return(p)
  }
)

#' Calculate cell trajectory in embedding space
#' @param embedding Original embedding matrix
#' @param delta_expr Expression changes
#' @param n_neighbors Number of neighbors for trajectory
#' @return Updated embedding with trajectory information
calculate_trajectory <- function(embedding, delta_expr, n_neighbors = 15) {
  if (!is.matrix(embedding) || !is.matrix(delta_expr)) {
    stop("Both embedding and delta_expr must be matrices")
  }

  nn <- RANN::nn2(embedding, k = n_neighbors)
  neighbors <- nn$nn.idx

  similarity_weights <- matrix(0, nrow = nrow(embedding), ncol = n_neighbors)
  for (i in seq_len(nrow(embedding))) {
    cell_neighbors <- neighbors[i, ]

    diff_vectors <- delta_expr[cell_neighbors, ] - matrix(
      delta_expr[i, ],
      nrow = length(cell_neighbors),
      ncol = ncol(delta_expr),
      byrow = TRUE
    )

    similarity_weights[i, ] <- exp(-rowSums(diff_vectors^2) / 0.1)
  }

  similarity_weights <- t(apply(similarity_weights, 1, function(x) x / sum(x)))

  trajectory_vectors <- matrix(0, nrow = nrow(embedding), ncol = ncol(embedding))
  for (i in seq_len(nrow(embedding))) {
    cell_neighbors <- neighbors[i, ]
    neighbor_coords <- embedding[cell_neighbors, ]
    current_coords <- matrix(embedding[i, ],
      nrow = length(cell_neighbors),
      ncol = ncol(embedding),
      byrow = TRUE
    )

    displacements <- neighbor_coords - current_coords
    trajectory_vectors[i, ] <- colSums(displacements * similarity_weights[i, ])
  }

  vector_lengths <- sqrt(rowSums(trajectory_vectors^2))
  valid_vectors <- vector_lengths > 0
  if (any(valid_vectors)) {
    trajectory_vectors[valid_vectors, ] <- sweep(
      trajectory_vectors[valid_vectors, ],
      1,
      vector_lengths[valid_vectors],
      "/"
    )
  }

  return(trajectory_vectors)
}
