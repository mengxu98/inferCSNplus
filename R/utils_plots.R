#' Plot goodness-of-fit metrics.
#'
#' @param object An object.
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param point_size Float indicating the point size.
#' @param ... Additional arguments.
#'
#' @return A ggplot2 object.
#'
#' @rdname plot_gof
#' @export
#' @method plot_gof CSNObject
setMethod(
  f = "plot_gof",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        point_size = 0.5,
                        ...) {
    if (is.null(celltypes)) {
      celltypes <- names(object@networks[[network]])
    } else {
      celltypes <- intersect(celltypes, names(object@networks[[network]]))
    }
    plot_list <- lapply(celltypes, function(celltype) {
      module_params <- NetworkModules(
        object,
        network = network
      )[[celltype]]@params

      if (length(module_params) == 0) {
        stop("No modules found, please run `find_modules()` first.")
      }

      metrics <- metrics(object, network = network)[[celltype]] |>
        dplyr::filter(r_squared <= 1, r_squared >= 0) |>
        dplyr::mutate(
          nice = r_squared > module_params$rsq_thresh & nvariables > module_params$nvar_thresh
        )

      p1 <- ggplot(
        metrics,
        aes(r_squared, nvariables, alpha = nice)
      ) +
        ggpointdensity::geom_pointdensity(size = point_size, shape = 16) +
        geom_hline(
          yintercept = module_params$nvar_thresh,
          size = 0.5,
          color = "darkgrey",
          linetype = "dashed"
        ) +
        geom_vline(
          xintercept = module_params$rsq_thresh,
          size = 0.5,
          color = "darkgrey",
          linetype = "dashed"
        ) +
        scale_color_gradientn(colors = rev(pals::brewer.reds(100))) +
        scale_alpha_discrete(range = c(0.5, 1)) +
        scale_y_continuous(
          trans = scales::pseudo_log_trans(base = 10),
          breaks = c(0, 1, 10, 100, 1000, 10000, 100000)
        ) +
        scale_x_continuous(breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        no_legend() +
        labs(
          x = expression("Explained variance" ~ (R**2)),
          y = "# variables in model"
        ) +
        theme(
          plot.margin = unit(c(0, 0, 0, 0), "line"),
          strip.text = element_blank()
        )

      p2 <- ggplot(metrics, aes(r_squared)) +
        geom_histogram(
          fill = "darkgray",
          bins = 20,
          color = "black",
          size = 0.2
        ) +
        theme_void() +
        no_legend()

      p3 <- ggplot(metrics, aes(nvariables)) +
        geom_histogram(
          fill = "darkgray",
          bins = 20,
          color = "black",
          size = 0.2
        ) +
        scale_x_continuous(
          trans = scales::pseudo_log_trans(base = 10),
          breaks = c(0, 1, 10, 100, 1000, 10000, 100000)
        ) +
        theme_void() +
        coord_flip() +
        no_legend()

      layout <- "
    AAAA#
    BBBBC
    BBBBC
    "
      p_out <- p2 + p1 + p3 + patchwork::plot_layout(design = layout) & no_margin()
    })
    names(plot_list) <- celltypes

    return(plot_list)
  }
)

#' Plot module metrics number of genes, number of peaks and number of TFs per gene.
#'
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#'
#' @return A ggplot2 object.
#'
#' @rdname plot_module_metrics
#' @export
#' @method plot_module_metrics CSNObject
setMethod(
  f = "plot_module_metrics",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    if (is.null(celltypes)) {
      celltypes <- names(object@networks[[network]])
    } else {
      celltypes <- intersect(
        celltypes,
        names(object@networks[[network]])
      )
    }

    plot_list <- purrr::map(
      celltypes, function(celltype) {
        modules <- NetworkModules(object, network = network)[[celltype]]@meta

        if (nrow(modules) == 0) {
          stop("No modules found, please run `find_modules()` first.")
        }

        plot_df <- modules |>
          dplyr::distinct(target, n_regions)

        p1 <- ggplot(plot_df, aes(1, n_regions)) +
          geom_violin(linewidth = 0.2, fill = "darkgrey", color = "black") +
          theme_bw() +
          no_x_text() +
          theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
          ) +
          geom_boxplot(
            width = 0.2,
            outlier.shape = NA,
            size = 0.2
          ) +
          labs(y = "# peaks") +
          ggtitle("# regions\nper target gene")


        plot_df <- modules |>
          dplyr::distinct(target, n_tfs)

        p2 <- ggplot(plot_df, aes(1, n_tfs)) +
          geom_violin(
            linewidth = 0.2,
            fill = "darkgrey",
            color = "black"
          ) +
          theme_bw() +
          no_x_text() +
          theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
          ) +
          geom_boxplot(
            width = 0.2,
            outlier.shape = NA,
            size = 0.2
          ) +
          labs(y = "# TFs") +
          ggtitle("# TFs\nper target gene")


        plot_df <- modules |>
          dplyr::distinct(tf, n_genes)

        p3 <- ggplot(plot_df, aes(1, n_genes)) +
          geom_violin(
            linewidth = 0.2,
            fill = "darkgrey",
            color = "black"
          ) +
          theme_bw() +
          no_x_text() +
          scale_y_continuous(
            trans = scales::pseudo_log_trans(base = 10)
          ) +
          theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
          ) +
          geom_boxplot(
            width = 0.2,
            outlier.shape = NA,
            size = 0.2
          ) +
          labs(y = expression("# genes")) +
          ggtitle("# target genes\nper TF")

        p_out <- p3 | p1 | p2 & no_margin()

        return(p_out)
      }
    )
    names(plot_list) <- celltypes

    return(plot_list)
  }
)

#' Compute UMAP embedding
#'
#' @param x data
#' @param n_pcs n_pcs
#' @param ... other parameters
#'
#' @return umap
get_umap <- function(
    x,
    n_pcs = 30,
    ...) {
  if (ncol(x) > 100) {
    pca_mat <- irlba::prcomp_irlba(x, n = n_pcs)$x
    rownames(pca_mat) <- rownames(x)
    x <- as.matrix(pca_mat)
  }
  umap_tbl <- uwot::umap(x, ...) %>%
    {
      colnames(.) <- c("UMAP_1", "UMAP_2")
      .
    } %>%
    tibble::as_tibble(rownames = "gene")

  return(umap_tbl)
}

#' Compute network graph embedding using UMAP.
#'
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param graph_name Name of the graph.
#' @param rna_assay Name of the RNA assay.
#' @param rna_layer Name of the RNA slot to use.
#' @param umap_method Method to compute edge weights for UMAP:
#' * \code{'weighted'} - Correlation weighted by GRN coefficient.
#' * \code{'corr'} - Only correlation.
#' * \code{'coef'} - Only GRN coefficient.
#' * \code{'none'} - Don't compute UMAP and create the graph directly
#' from modules.
#' @param features Features to use to create the graph. If \code{NULL}
#' uses all features in the network.
#' @param seed Random seed for UMAP computation
#' @param ... Additional arguments for \code{\link[uwot]{umap}}.
#' @param verbose Print messages.
#'
#' @return A CSNObject object.
#'
#' @rdname get_network_graph
#' @export
#' @method get_network_graph CSNObject
setMethod(
  f = "get_network_graph",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        graph_name = "module_graph",
                        rna_assay = "RNA",
                        rna_layer = "data",
                        umap_method = c("weighted", "corr", "coef", "none"),
                        features = NULL,
                        seed = 1,
                        verbose = TRUE,
                        ...) {
    umap_method <- match.arg(umap_method)
    modules_list <- NetworkModules(object, network = network)

    if (is.null(celltypes)) {
      celltypes <- names(modules_list)
    } else {
      celltypes <- intersect(celltypes, names(modules_list))
    }
    gene_graph_list <- purrr::map(
      celltypes, function(celltype) {
        modules <- modules_list[[celltype]]
        if (length(modules@params) == 0) {
          stop("No modules found, please run `find_modules()` first.")
        }

        if (is.null(features)) {
          features <- get_attribute(
            object,
            attribute = "targets",
            active_network = network,
            celltypes = celltype
          )
        }

        if (umap_method == "weighted") {
          rna_expr <- t(
            LayerData(
              object,
              assay = rna_assay,
              layer = rna_layer
            )
          )
          features <- intersect(features, colnames(rna_expr))

          log_message("Computing gene-gene correlation", verbose = verbose)
          rna_expr <- rna_expr[, features]
          gene_cor <- sparse_cor(rna_expr)
          gene_cor_df <- gene_cor |>
            tibble::as_tibble(rownames = "source") |>
            tidyr::pivot_longer(
              !source,
              names_to = "target",
              values_to = "corr"
            )

          modules_use <- modules@meta |>
            dplyr::filter(
              target %in% colnames(rna_expr),
              tf %in% colnames(rna_expr)
            )

          gene_net <- modules_use %>%
            dplyr::select(tf, target, tidyselect::everything()) %>%
            dplyr::group_by(target) %>%
            left_join(
              gene_cor_df,
              by = c("tf" = "source", "target")
            ) %>%
            {
              .$corr[is.na(.$corr)] <- 0
              .
            }

          reg_mat <- gene_net %>%
            dplyr::select(target, tf, coefficient) %>%
            tidyr::pivot_wider(
              names_from = tf,
              values_from = coefficient,
              values_fill = 0
            ) %>%
            tibble::column_to_rownames("target") %>%
            as.matrix()

          log_message("Computing weighted regulatory factor", verbose = verbose)
          reg_factor_mat <- abs(reg_mat) + 1
          coex_mat <- gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)
        } else if (umap_method == "corr") {
          net_features <- get_attribute(
            object,
            celltypes = NULL,
            attribute = "targets",
            active_network = network
          )
          rna_expr <- t(LayerData(object, assay = rna_assay, layer = rna_layer))

          if (!is.null(features)) {
            features <- intersect(
              intersect(features, colnames(rna_expr)),
              net_features
            )
          } else {
            features <- net_features
          }

          log_message("Computing gene-gene correlation", verbose = verbose)
          rna_expr <- rna_expr[, features]
          coex_mat <- sparse_cor(rna_expr)
          gene_cor_df <- coex_mat %>%
            tibble::as_tibble(rownames = "source") %>%
            tidyr::pivot_longer(
              !source,
              names_to = "target",
              values_to = "corr"
            )

          # Get adjacency df and matrix
          modules_use <- modules@meta %>%
            dplyr::filter(target %in% features, tf %in% features)

          gene_net <- modules_use %>%
            dplyr::select(tf, target, tidyselect::everything()) %>%
            dplyr::group_by(target) %>%
            left_join(gene_cor_df, by = c("tf" = "source", "target")) %>%
            {
              .$corr[is.na(.$corr)] <- 0
              .
            }
        } else if (umap_method == "coef") {
          modules_use <- modules@meta %>%
            dplyr::filter(target %in% features, tf %in% features)

          gene_net <- modules_use %>%
            dplyr::select(tf, target, tidyselect::everything()) %>%
            dplyr::group_by(target)

          coex_mat <- gene_net %>%
            dplyr::select(target, tf, coefficient) %>%
            tidyr::pivot_wider(
              names_from = tf,
              values_from = coefficient,
              values_fill = 0
            ) %>%
            tibble::column_to_rownames("target") %>%
            as.matrix()
        } else if (umap_method == "none" | is.null(umap_method)) {
          modules_use <- modules@meta %>%
            dplyr::filter(target %in% features, tf %in% features)

          gene_net <- modules_use %>%
            dplyr::select(tf, target, tidyselect::everything()) %>%
            dplyr::group_by(target)

          log_message("Getting network graph", verbose = verbose)
          gene_graph <- tidygraph::as_tbl_graph(gene_net) %>%
            tidygraph::activate(edges) %>%
            dplyr::mutate(
              from_node = tidygraph::.N()$name[from],
              to_node = tidygraph::.N()$name[to]
            ) %>%
            dplyr::mutate(dir = sign(coefficient)) %>%
            tidygraph::activate(nodes) %>%
            dplyr::mutate(centrality = tidygraph::centrality_pagerank())

          object@csn@networks[[network]]@graphs[[graph_name]] <- gene_graph
          return(object)
        }

        log_message("Computing UMAP embedding", verbose = verbose)
        set.seed(seed)
        coex_umap <- get_umap(
          as_matrix(coex_mat),
          ...
        )

        log_message("Getting network graph", verbose = verbose)
        gene_graph <- tidygraph::as_tbl_graph(gene_net) |>
          tidygraph::activate(edges) |>
          dplyr::mutate(
            from_node = tidygraph::.N()$name[from],
            to_node = tidygraph::.N()$name[to]
          ) |>
          dplyr::mutate(dir = sign(coefficient)) |>
          tidygraph::activate(nodes) |>
          dplyr::mutate(centrality = tidygraph::centrality_pagerank()) |>
          dplyr::inner_join(coex_umap, by = c("name" = "gene"))
      }
    )
    names(gene_graph_list) <- celltypes
    for (celltype in celltypes) {
      object@networks[[network]][[celltype]]@graphs[[graph_name]] <- gene_graph_list[[celltype]]
    }

    return(object)
  }
)

#' Plot network graph.
#'
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param graph Name of the graph.
#' @param layout Layout for the graph. Can be 'umap' or any force-directed layout
#' implemented in \code{\link[ggraph]{ggraph}}
#' @param edge_width Edge width.
#' @param edge_color Edge color.
#' @param node_color Node color or color gradient.
#' @param node_size Node size range.
#' @param text_size Font size for labels.
#' @param color_nodes Logical, Whether to color nodes by centrality.
#' @param label_nodes Logical, Whether to label nodes with gene name.
#' @param color_edges Logical, Whether to color edges by direction.
#'
#' @return A CSNObject object.
#'
#' @rdname plot_network_graph
#' @export
#' @method plot_network_graph CSNObject
setMethod(
  f = "plot_network_graph",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        graph = "module_graph",
                        layout = "umap",
                        edge_width = 0.2,
                        edge_color = c("-1" = "darkgrey", "1" = "orange"),
                        node_color = pals::magma(100),
                        node_size = c(1, 5),
                        text_size = 10,
                        color_nodes = TRUE,
                        label_nodes = TRUE,
                        color_edges = TRUE,
                        ...) {
    gene_graph_list <- NetworkGraph(
      object,
      network = network,
      graph = graph
    )
    if (is.null(celltypes)) {
      celltypes <- names(gene_graph_list)
    } else {
      celltypes <- intersect(celltypes, names(gene_graph_list))
    }
    plot_list <- purrr::map(
      celltypes,
      function(celltype) {
        gene_graph <- gene_graph_list[[celltype]]
        has_umap <- "UMAP_1" %in%
          colnames(tibble::as_tibble(tidygraph::activate(gene_graph, "nodes")))
        if (layout == "umap" && !has_umap) {
          stop("No UMAP coordinates found, please run `get_network_graph()` first.")
        }

        if (layout == "umap") {
          p <- ggraph::ggraph(gene_graph, x = UMAP_1, y = UMAP_2)
        } else {
          p <- ggraph::ggraph(gene_graph, layout = layout)
        }

        if (color_edges) {
          p <- p + ggraph::geom_edge_diagonal(
            aes(color = factor(dir)),
            width = edge_width
          ) +
            ggraph::scale_edge_color_manual(values = edge_color)
        } else {
          p <- p + ggraph::geom_edge_diagonal(
            width = edge_width,
            color = edge_color[1]
          )
        }

        if (color_nodes) {
          p <- p +
            ggraph::geom_node_point(
              aes(fill = centrality, size = centrality),
              color = "darkgrey",
              shape = 21
            ) +
            scale_fill_gradientn(colors = node_color)
        } else {
          p <- p + ggraph::geom_node_point(
            color = "darkgrey",
            shape = 21,
            fill = "lightgrey",
            size = node_size[1],
            stroke = 0.5
          )
        }

        if (label_nodes) {
          p <- p + ggraph::geom_node_text(
            aes(label = name),
            repel = TRUE,
            size = text_size / ggplot2::.pt,
            max.overlaps = 99999
          )
        }
        p <- p + scale_size_continuous(range = node_size) +
          theme_void() + no_legend()
      }
    )
    names(plot_list) <- celltypes

    return(plot_list)
  }
)

#' Get sub-network centered around one TF.
#'
#' @param tfs The transcription factors to center around.
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param graph Name of the graph.
#' @param features Features to use. If \code{NULL} uses all features in the graph.
#' @param order Integer indicating the maximal order of the graph.
#' @param keep_all_edges Logical, whether to maintain all edges to each leaf
#' or prune to the strongest overall connection.
#' @param verbose Logical. Whether to print messages.
#' @param cores Logical. Whether to parallelize the computation with \code{\link[foreach]{foreach}}.
#'
#' @return A CSNObject object.
#'
#' @rdname get_tf_network
#' @export
#' @method get_tf_network CSNObject
setMethod(
  f = "get_tf_network",
  signature = "CSNObject",
  definition = function(object,
                        celltypes = NULL,
                        tfs = NULL,
                        features = NULL,
                        network = DefaultNetwork(object),
                        graph = "module_graph",
                        order = 3,
                        keep_all_edges = FALSE,
                        verbose = TRUE,
                        cores = 1,
                        ...) {
    gene_graph_list <- NetworkGraph(object, network = network, graph = graph)
    if (is.null(celltypes)) {
      celltypes <- names(gene_graph_list)
    } else {
      celltypes <- intersect(celltypes, names(gene_graph_list))
    }

    if (is.null(tfs)) {
      tfs <- get_attribute(
        object,
        attribute = "regulators",
        active_network = network
      )
    }

    tfs_celltype_graph_list <- purrr::map(
      celltypes,
      function(celltype) {
        tfs_graph_list <- purrr::map(
          tfs,
          function(tf) {
            gene_graph <- gene_graph_list[[celltype]]

            gene_graph_nodes <- gene_graph |>
              tidygraph::activate(nodes) |>
              tibble::as_tibble()

            if (is.null(features)) {
              features <- get_attribute(
                object,
                celltypes = celltype,
                attribute = "targets",
                active_network = network
              )
            }

            features <- intersect(features, gene_graph_nodes$name)

            log_message("Getting shortest paths from TF", verbose = verbose)
            spaths <- igraph::all_shortest_paths(
              gene_graph,
              tf,
              features,
              mode = "out"
            )$res

            if (length(spaths) < 3) {
              log_message(
                "Selected TF has fewer than 3 targets.",
                message_type = "warning",
                verbose = verbose
              )
              return(NULL)
            }

            spath_list <- parallelize_fun(spaths, function(p) {
              edg <- names(p)
              edg_graph <- gene_graph |>
                dplyr::filter(name %in% edg) |>
                tidygraph::convert(
                  to_shortest_path,
                  from = which(tidygraph::.N()$name == edg[1]),
                  to = which(tidygraph::.N()$name == edg[length(edg)])
                ) %E>%
                dplyr::mutate(
                  from_node = tidygraph::.N()$name[from],
                  to_node = tidygraph::.N()$name[to]
                ) |>
                tibble::as_tibble()

              edg_dir <- edg_graph |>
                dplyr::pull(coefficient) |>
                sign() |>
                prod()
              edg_est <- edg_graph |>
                dplyr::pull(coefficient) |>
                mean()
              path_df <- tibble::tibble(
                start_node = edg[1],
                end_node = edg[length(edg)],
                dir = edg_dir,
                path = paste(edg, collapse = ";"),
                path_regions = paste(edg_graph$regions, collapse = ";"),
                order = length(edg) - 1,
                mean_estimate = edg_est
              )

              if ("padj" %in% colnames(edg_graph)) {
                path_df$mean_log_padj <- edg_graph |>
                  dplyr::pull(padj) %>%
                  {
                    -log10(.)
                  } |>
                  mean()
              }

              return(
                list(
                  path = path_df,
                  graph = dplyr::mutate(
                    edg_graph,
                    path = paste(edg, collapse = ";"),
                    end_node = edg[length(edg)],
                    comb_dir = edg_dir
                  )
                )
              )
            }, cores = cores, verbose = verbose)

            log_message("Pruning graph", verbose = verbose)
            spath_dir <- purrr::map_dfr(spath_list, function(x) x$path) |>
              dplyr::mutate(
                path_genes = stringr::str_split(path, ";"),
                path_regions = stringr::str_split(path_regions, ";")
              )
            spath_graph <- purrr::map_dfr(spath_list, function(x) x$graph)

            csn_pruned <- spath_dir |>
              dplyr::select(
                start_node,
                end_node,
                tidyselect::everything()
              ) |>
              dplyr::group_by(end_node) |>
              dplyr::filter(order <= order)

            if (!keep_all_edges) {
              if ("padj" %in% colnames(csn_pruned)) {
                csn_pruned <- dplyr::filter(
                  csn_pruned,
                  order == 1 | mean_padj == max(mean_padj)
                )
              } else {
                csn_pruned <- dplyr::filter(
                  csn_pruned,
                  order == 1 | mean_estimate == max(mean_estimate)
                )
              }
            }

            spath_graph_pruned <- spath_graph |>
              dplyr::filter(path %in% csn_pruned$path) |>
              dplyr::select(from_node, to_node, end_node, comb_dir) |>
              dplyr::distinct()

            csn_graph_pruned <- gene_graph %E>%
              dplyr::mutate(
                from_node = tidygraph::.N()$name[from],
                to_node = tidygraph::.N()$name[to]
              ) |>
              tibble::as_tibble() |>
              dplyr::distinct()
            csn_graph_pruned <- suppressMessages(
              inner_join(csn_graph_pruned, spath_graph_pruned)
            ) |>
              dplyr::select(
                from_node,
                to_node,
                tidyselect::everything(),
                -from,
                -to
              ) |>
              dplyr::arrange(comb_dir) |>
              tidygraph::as_tbl_graph()

            return(csn_graph_pruned)
          }
        )
        names(tfs_graph_list) <- tfs
        tfs_graph_list <- tfs_graph_list[!purrr::map_lgl(tfs_graph_list, is.null)]

        return(tfs_graph_list)
      }
    )
    names(tfs_celltype_graph_list) <- celltypes
    for (celltype in celltypes) {
      object@networks[[network]][[celltype]]@graphs$tf_graphs <- tfs_celltype_graph_list[[celltype]]
    }
    return(object)
  }
)


#' Plot sub-network centered around one TF.
#'
#' @param tfs The transcription factor to center around.
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param graph Name of the graph.
#' @param circular Logical. Layout tree in circular layout.
#' @param edge_width Edge width.
#' @param edge_color Edge color.
#' @param node_size Node size.
#' @param text_size Font size for labels.
#' @param label_nodes String, indicating what to label.
#' * \code{'tfs'} - Label all TFs.
#' * \code{'all'} - Label all genes.
#' * \code{'none'} - Label nothing (except the root TF).
#' @param color_edges Logical, whether to color edges by direction.
#'
#' @return A CSNObject object.
#'
#' @rdname plot_tf_network
#' @export
#' @method plot_tf_network CSNObject
setMethod(
  f = "plot_tf_network",
  signature = "CSNObject",
  definition = function(object,
                        tfs = NULL,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        graph = "module_graph",
                        circular = TRUE,
                        edge_width = 0.2,
                        edge_color = c("-1" = "darkgrey", "1" = "orange"),
                        node_size = 3,
                        text_size = 10,
                        label_nodes = c("tfs", "all", "none"),
                        color_edges = TRUE,
                        ...) {
    label_nodes <- match.arg(label_nodes)
    if (is.null(celltypes)) {
      celltypes <- get_attribute(
        object,
        attribute = "celltypes",
        active_network = network
      )
    } else {
      celltypes <- intersect(celltypes, get_attribute(
        object,
        attribute = "celltypes",
        active_network = network
      ))
    }

    gene_graph_list <- NetworkGraph(
      object,
      network = network,
      graph = "tf_graphs",
      celltypes = celltypes
    )
    names(gene_graph_list) <- celltypes

    plot_list_celltype <- purrr::map(
      celltypes,
      function(celltype) {
        if (is.null(tfs)) {
          tfs <- get_attribute(
            object,
            attribute = "regulators",
            active_network = network,
            celltypes = celltype
          )
        } else {
          tfs <- intersect(tfs, get_attribute(
            object,
            attribute = "regulators",
            active_network = network,
            celltypes = celltype
          ))
        }
        gene_graph_celltype <- gene_graph_list[[celltype]]
        tfs <- intersect(tfs, names(gene_graph_celltype))

        plot_list <- purrr::map(tfs, function(tf) {
          gene_graph <- gene_graph_celltype[[tf]]
          if (is.null(gene_graph)) {
            log_message(
              paste0("No graph found for TF ", tf, " in celltype ", celltype),
              message_type = "warning"
            )
            return(NULL)
          }

          p <- ggraph::ggraph(
            gene_graph,
            layout = "tree",
            circular = circular
          )

          if (color_edges) {
            p <- p +
              ggraph::geom_edge_diagonal(
                aes(color = factor(dir)),
                width = edge_width
              ) +
              ggraph::scale_edge_color_manual(
                values = edge_color
              )
          } else {
            p <- p +
              ggraph::geom_edge_diagonal(
                width = edge_width, color = edge_color[1]
              )
          }

          p <- p +
            ggraph::geom_node_point(
              color = "darkgrey",
              shape = 21,
              fill = "lightgrey",
              size = node_size,
              stroke = 0.5
            )

          net_tfs <- colnames(NetworkTFs(object))
          if (label_nodes == "tfs") {
            p <- p +
              ggraph::geom_node_label(
                aes(label = name, filter = name %in% net_tfs),
                size = text_size / ggplot2::.pt,
                label.padding = unit(0.1, "line")
              )
          } else if (label_nodes == "all") {
            p <- p +
              ggraph::geom_node_label(
                aes(label = name),
                size = text_size / ggplot2::.pt,
                label.padding = unit(0.1, "line")
              )
          } else {
            p <- p +
              ggraph::geom_node_label(
                aes(label = name, filter = name == tf),
                size = text_size / ggplot2::.pt,
                label.padding = unit(0.1, "line")
              )
          }
          p <- p + scale_size_continuous(range = node_size) +
            theme_void() + no_legend()
        })
        names(plot_list) <- tfs
        plot_list <- plot_list[!purrr::map_lgl(plot_list, is.null)]

        return(plot_list)
      }
    )
    names(plot_list_celltype) <- celltypes

    return(plot_list_celltype)
  }
)
