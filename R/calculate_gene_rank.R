#' @include setClass.R
#' @include setGenerics.R

#' @title Calculate gene rank
#' @param object Network object
#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @param ... Additional arguments
#' @return Data frame with gene ranks
#' @export
#' @rdname calculate_gene_rank
setGeneric(
  name = "calculate_gene_rank",
  signature = c("object"),
  def = function(object,
                 regulators = NULL,
                 targets = NULL,
                 directed = FALSE,
                 method = c("page_rank", "degree_distribution"),
                 ...) {
    standardGeneric("calculate_gene_rank")
  }
)

#' @rdname calculate_gene_rank
setMethod(
  "calculate_gene_rank",
  signature(object = "Network"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE,
           method = c("page_rank", "degree_distribution"),
           ...) {
    method <- match.arg(method)
    if (method == "page_rank") {
      return(calculate_page_rank(object, regulators, targets, directed, ...))
    } else {
      return(calculate_degree_distribution(object, regulators, targets, directed, ...))
    }
  }
)

#' @rdname calculate_gene_rank
setMethod(
  "calculate_gene_rank",
  signature(object = "data.frame"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE,
           method = c("page_rank", "degree_distribution"),
           ...) {
    method <- match.arg(method)
    if (method == "page_rank") {
      return(calculate_page_rank(object, regulators, targets, directed, ...))
    } else {
      return(calculate_degree_distribution(object, regulators, targets, directed, ...))
    }
  }
)

#' @title Calculate PageRank
#' @param object Network object
#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @return Data frame with PageRank scores
#' @export
#' @rdname calculate_page_rank
setGeneric(
  name = "calculate_page_rank",
  signature = c("object"),
  def = function(object,
                 regulators = NULL,
                 targets = NULL,
                 directed = FALSE) {
    standardGeneric("calculate_page_rank")
  }
)

#' @rdname calculate_page_rank
setMethod(
  "calculate_page_rank",
  signature(object = "Network"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- as.data.frame(object@network)
    .calculate_page_rank(network_table, directed)
  }
)

#' @rdname calculate_page_rank
setMethod(
  "calculate_page_rank",
  signature(object = "data.frame"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- network_format(
      object,
      regulators,
      targets,
      abs_weight = FALSE
    )
    .calculate_page_rank(network_table, directed)
  }
)

.calculate_page_rank <- function(network_table, directed) {
  network <- igraph::graph_from_data_frame(
    network_table,
    directed = directed
  )
  page_rank_res <- data.frame(
    igraph::page_rank(network, directed = directed)$vector
  )
  colnames(page_rank_res) <- c("page_rank")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "page_rank")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$page_rank,
    decreasing = TRUE
  ), ]
  page_rank_res$regulator <- ifelse(
    page_rank_res$gene %in% unique(network_table$regulator),
    "TRUE", "FALSE"
  )
  rownames(page_rank_res) <- NULL

  return(page_rank_res)
}

#' @rdname calculate_gene_rank
setGeneric(
  name = "calculate_degree_distribution",
  signature = c("object"),
  def = function(object,
                 regulators = NULL,
                 targets = NULL,
                 directed = TRUE) {
    standardGeneric("calculate_degree_distribution")
  }
)

#' @rdname calculate_gene_rank
setMethod(
  "calculate_degree_distribution",
  signature(object = "Network"),
  function(object, regulators = NULL, targets = NULL, directed = TRUE) {
    network_table <- as.data.frame(object@network)
    if (!all(c("regulator", "target", "weight") %in% colnames(network_table))) {
      stop("Network slot must contain 'regulator', 'target', and 'weight' columns")
    }
    .calculate_degree_distribution_internal(network_table, directed)
  }
)

#' @rdname calculate_gene_rank
setMethod(
  "calculate_degree_distribution",
  signature(object = "data.frame"),
  function(object, regulators = NULL, targets = NULL, directed = TRUE) {
    network_table <- network_format(
      object,
      regulators,
      targets,
      abs_weight = FALSE
    )
    .calculate_degree_distribution_internal(network_table, directed)
  }
)

.calculate_power_law_fit <- function(deg) {
  degree_freq <- table(deg)
  k <- as.numeric(names(degree_freq))
  pk <- as.numeric(degree_freq) / sum(degree_freq)

  if (length(unique(k)) > 1) {
    log_k <- log(k)
    log_pk <- log(pk)
    fit <- lm(log_pk ~ log_k)
    return(summary(fit)$r.squared)
  }
  return(0)
}

.calculate_degree_distribution_internal <- function(
    network_table, directed) {
  network <- igraph::graph_from_data_frame(
    network_table,
    directed = directed
  )

  total_degrees <- igraph::degree(network, mode = "total")

  result <- data.frame(
    gene = names(total_degrees),
    degree = as.numeric(total_degrees)
  )

  degree_freq <- table(total_degrees)
  result$P_k <- degree_freq[match(result$degree, names(degree_freq))] / sum(degree_freq)

  if (directed) {
    in_degrees <- igraph::degree(network, mode = "in")
    out_degrees <- igraph::degree(network, mode = "out")

    in_power_law_score <- .calculate_power_law_fit(in_degrees)
    out_power_law_score <- .calculate_power_law_fit(out_degrees)

    result$in_degree <- in_degrees[result$gene]
    result$out_degree <- out_degrees[result$gene]
    result$in_power_law_fit <- in_power_law_score
    result$out_power_law_fit <- out_power_law_score
    result$rank_value <- (result$in_degree * in_power_law_score +
      result$out_degree * out_power_law_score) / 2
  } else {
    power_law_score <- .calculate_power_law_fit(total_degrees)
    result$power_law_fit <- power_law_score
    result$rank_value <- result$degree * power_law_score
  }

  result <- result[order(result$rank_value, decreasing = TRUE), ]
  result$regulator <- result$gene %in% unique(network_table$regulator)
  rownames(result) <- NULL

  return(result)
}

#' @title Plot gene ranks and network properties
#' @param object Network object
#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @param weight_cutoff Numeric, threshold for edge weight filtering
#' @param compare_random Logical, whether to compare with randomized network
#' @return Combined ggplot object
#' @export
#' @rdname plot_gene_rank
setGeneric(
  name = "plot_gene_rank",
  signature = c("object"),
  def = function(object,
                method = c("page_rank", "degree_distribution"),
                weight_cutoff = 0.1,
                compare_random = TRUE) {
    standardGeneric("plot_gene_rank")
  }
)

#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "Network"),
  function(object,
           method = c("page_rank", "degree_distribution"),
           weight_cutoff = 0.1,
           compare_random = TRUE) {
    method <- match.arg(method)
    network_table <- as.data.frame(object@network)
    network_table <- network_table[abs(network_table$weight) >= weight_cutoff, ]

    g_orig <- igraph::graph_from_data_frame(network_table, directed = TRUE)
    degrees_orig <- igraph::degree(g_orig, mode = "total")

    degree_freq <- table(degrees_orig)
    df_orig <- data.frame(
      k = as.numeric(names(degree_freq)),
      P_k = as.numeric(degree_freq) / sum(degree_freq)
    )

    dist_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey95", size = 0.2),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 10),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        aspect.ratio = 1
      )

    p1 <- ggplot(df_orig, aes(x = log(k), y = log(P_k))) +
      geom_point(size = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "grey50", size = 0.3, formula = y ~ x) +
      dist_theme +
      labs(x = "log k", y = "log P(k)", title = "Degree distribution", tag = "a")

    if (nrow(df_orig) > 1) {
      model <- lm(log(P_k) ~ log(k), data = df_orig)
      r2 <- summary(model)$r.squared
      p1 <- p1 + annotate("text",
        x = max(log(df_orig$k)) - 0.1,
        y = max(log(df_orig$P_k)) - 0.1,
        label = sprintf("R² = %.3f", r2),
        size = 3,
        color = "steelblue",
        hjust = 1
      )
    }

    if (method == "page_rank") {
      ranks <- calculate_page_rank(object)
      rank_col <- "page_rank"
    } else {
      ranks <- calculate_degree_distribution(object)
      rank_col <- "rank_value"
    }

    centrality_df <- ranks %>%
      arrange(desc(.data[[rank_col]])) %>%
      head(30) %>%
      mutate(centrality = .data[[rank_col]] / max(.data[[rank_col]]))

    cent_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey95", size = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3)
      )

    p_cent <- ggplot(
      centrality_df,
      aes(x = centrality, y = reorder(gene, centrality))
    ) +
      geom_point(color = "steelblue", size = 2) +
      cent_theme +
      labs(
        x = if (method == "page_rank") "PageRank centrality" else "Degree centrality",
        y = NULL,
        tag = "c"
      )

    if (compare_random) {
      n_edges <- igraph::ecount(g_orig)
      all_nodes <- unique(c(network_table$regulator, network_table$target))
      n_nodes <- length(all_nodes)

      random_edges <- data.frame(
        regulator = sample(all_nodes, n_edges, replace = TRUE),
        target = sample(all_nodes, n_edges, replace = TRUE),
        weight = network_table$weight
      )

      g_random <- igraph::graph_from_data_frame(random_edges, directed = TRUE)
      degrees_random <- igraph::degree(g_random, mode = "total")

      degree_freq_random <- table(degrees_random)
      df_random <- data.frame(
        k = as.numeric(names(degree_freq_random)),
        P_k = as.numeric(degree_freq_random) / sum(degree_freq_random)
      )

      p2 <- ggplot(df_random, aes(x = log(k), y = log(P_k))) +
        geom_point(size = 1) +
        geom_smooth(method = "lm", se = FALSE, color = "grey50", size = 0.3, formula = y ~ x) +
        dist_theme +
        labs(
          x = "log k", y = "log P(k)",
          title = "Degree distribution\nof randomized network",
          tag = "b"
        )

      if (nrow(df_random) > 1) {
        model_random <- lm(log(P_k) ~ log(k), data = df_random)
        r2_random <- summary(model_random)$r.squared
        p2 <- p2 + annotate("text",
          x = max(log(df_random$k)) - 0.1,
          y = max(log(df_random$P_k)) - 0.1,
          label = sprintf("R² = %.3f", r2_random),
          size = 3,
          color = "steelblue",
          hjust = 1
        )
      }
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (compare_random) {
        dist_plots <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1))
        return(dist_plots | p_cent + patchwork::plot_layout(widths = c(1, 1.5)))
      }
      return(p1 | p_cent + patchwork::plot_layout(widths = c(1, 1.5)))
    }

    return(
      list(
        distribution = if (compare_random) p1 / p2 else p1,
        centrality = p_cent
      )
    )
  }
)

#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "data.frame"),
  function(object,
           method = c("page_rank", "degree_distribution"),
           weight_cutoff = 0.1,
           compare_random = TRUE) {
    network_obj <- new("Network", network = object)
    plot_gene_rank(network_obj, method, weight_cutoff, compare_random)
  }
)
