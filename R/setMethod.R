#' @title Get GRN inference parameters
#' @rdname Params
#' @export
setMethod(
  f = "Params",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@params)
  }
)

#' @include setClass.R
#' @include setGenerics.R

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname GetNetwork
#' @export
setMethod(
  f = "GetNetwork",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    if (!is.null(celltypes)) {
      return(
        lapply(
          celltypes,
          function(x) {
            object@networks[[network]][[x]]
          }
        )
      )
    }
    return(object@networks[[network]])
  }
)

#' @rdname NetworkTFs
#' @export
setMethod(
  f = "NetworkTFs",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@regions@motifs2tfs)
  }
)

#' @rdname NetworkRegions
#' @export
setMethod(
  f = "NetworkRegions",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@regions)
  }
)

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkModules
#' @export
setMethod(
  f = "NetworkModules",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkModules(net)
        )
      )
    }
    return(
      lapply(
        networks,
        function(net) NetworkModules(net)
      )
    )
  }
)

#' @rdname NetworkModules
#' @export
setMethod(
  f = "NetworkModules",
  signature = "Network",
  definition = function(object, ...) {
    return(object@modules)
  }
)

#' @rdname NetworkParams
#' @export
setMethod(
  f = "NetworkParams",
  signature = "Network",
  definition = function(object, ...) {
    return(object@params)
  }
)

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkParams
#' @export
setMethod(
  f = "NetworkParams",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkParams(net)
        )
      )
    }
    return(
      lapply(
        networks, function(net) NetworkParams(net)
      )
    )
  }
)

#' @param network network
#' @param graph graph
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkGraph
#' @export
setMethod(
  f = "NetworkGraph",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        graph = "module_graph",
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkGraph(net, graph = graph)
        )
      )
    }
    return(
      lapply(
        networks,
        function(net) NetworkGraph(net, graph = graph)
      )
    )
  }
)

#' @rdname NetworkGraph
#' @export
setMethod(
  f = "NetworkGraph",
  signature = "Network",
  definition = function(object,
                        graph = "module_graph",
                        ...) {
    if (!graph %in% names(object@graphs)) {
      stop(
        paste0(
          "The requested graph '", graph, "' does not exist. ",
          "Try (re-)running `get_network_graph()`."
        )
      )
    }
    return(object@graphs[[graph]])
  }
)

#' @title Get active network
#'
#' @rdname DefaultNetwork
#' @method DefaultNetwork CSNObject
#' @export
setMethod(
  f = "DefaultNetwork",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@active_network)
  }
)

#' @param group_name group_name
#' @param assay assay
#' @param verbose verbose
#'
#' @rdname GetAssaySummary
#' @export
setMethod(
  f = "GetAssaySummary",
  signature = "Seurat",
  definition = function(object,
                        group_name,
                        assay = NULL,
                        verbose = TRUE,
                        ...) {
    if (is.null(assay)) {
      assay <- object@active.assay
    }
    smry <- Seurat::Misc(
      object[[assay]]
    )$summary[[group_name]]
    if (is.null(smry)) {
      log_message(
        "summary of '", group_name, "' does not yet exist",
        verbose = verbose,
        message_type = "warning"
      )
      log_message(
        "Summarizing information for '", group_name, "'",
        verbose = verbose
      )
      object <- aggregate_assay(
        object,
        assay = assay,
        group_name = group_name
      )
      smry <- GetAssaySummary(
        object,
        assay = assay,
        group_name = group_name,
        verbose = verbose
      )
    }
    return(smry)
  }
)

#' @title Get summary of seurat assay from CSNObject
#'
#' @param group_name group_name
#' @param assay assay
#' @param verbose verbose
#'
#' @rdname GetAssaySummary
#' @export
setMethod(
  f = "GetAssaySummary",
  signature = "CSNObject",
  definition = function(object,
                        group_name,
                        assay = NULL,
                        verbose = TRUE,
                        ...) {
    return(
      GetAssaySummary(
        object@data,
        group_name,
        assay = assay,
        verbose = verbose
      )
    )
  }
)

#' @title Get Seurat assay from CSNObject
#'
#' @param assay assay
#'
#' @rdname GetAssay
#' @method GetAssay CSNObject
#' @export
setMethod(
  f = "GetAssay",
  signature = "CSNObject",
  definition = function(object,
                        assay = NULL,
                        ...) {
    return(
      Seurat::GetAssay(
        object@data,
        assay = assay
      )
    )
  }
)

#' @rdname LayerData
#' @export
setMethod(
  f = "LayerData",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(SeuratObject::LayerData(object@data, ...))
  }
)

#' @rdname VariableFeatures
#' @export
setMethod(
  f = "VariableFeatures",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(Seurat::VariableFeatures(object@data, ...))
  }
)

#' @title Get fitted coefficients
#'
#' @param object Network object
#' @param ... Other parameters
#'
#' @method coef Network
#'
#' @return Return the fitted coefficients
#' @export
coef.Network <- function(object, ...) {
  return(object@coefficients)
}

#' @title Get fitted coefficients
#'
#' @param object CSNObject object
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' If \code{NULL}, all celltypes are plotted.
#' @param ... Other parameters
#'
#' @method coef CSNObject
#'
#' @return Return the fitted coefficients
#' @export
coef.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    celltypes = NULL,
    ...) {
  networks <- GetNetwork(
    object,
    network = network,
    celltypes = celltypes
  )
  if (is.null(celltypes)) {
    return(
      lapply(
        networks, function(net) {
          net@coefficients
        }
      )
    )
  }
  return(networks@coefficients)
}

#' @title Get goodness-of-fit info
#'
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
# metrics methods
#' @rdname metrics
#' @export
setMethod(
  f = "metrics",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    celltypes_all <- get_attribute(
      object,
      attribute = "celltypes"
    )
    celltypes <- intersect(
      celltypes %ss% celltypes_all,
      celltypes_all
    )
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )

    res <- lapply(
      networks,
      function(net) {
        net@metrics
      }
    )
    names(res) <- celltypes

    return(res)
  }
)

#' @rdname metrics
#' @export
setMethod(
  f = "metrics",
  signature = "Network",
  definition = function(object, celltypes = NULL, ...) {
    return(object@metrics)
  }
)

#' @title Print Network objects
#'
#' @param x x
#' @param ... other parameters
#'
#' @rdname print
#' @export
#' @method print Network
print.Network <- function(x, ...) {
  coeffs <- methods::slot(x, "coefficients")
  modules_obj <- methods::slot(x, "modules")
  modules_meta <- methods::slot(modules_obj, "meta")

  if (nrow(modules_meta) == 0) {
    n_genes <- length(unique(coeffs$target))
    n_tfs <- length(unique(coeffs$tf))
  } else {
    n_genes <- length(unique(modules_meta$target))
    n_tfs <- length(unique(modules_meta$tf))
  }
  cat(paste0(
    "A Network object\n", "with ", n_tfs, " TFs and ",
    n_genes, " target genes"
  ))
}

setMethod(
  "show",
  signature = "Network",
  function(object) {
    print(object)
  }
)

#' @title Print Modules objects
#'
#' @rdname print
#' @export
#' @method print Modules
print.Modules <- function(x, ...) {
  n_mods <- length(x@features$genes_pos)
  cat(paste0(
    "An Modules object with ", n_mods, " TF modules"
  ))
}

setMethod(
  "show",
  signature = "Modules",
  function(object) {
    print(object)
  }
)

#' @title Print Regions objects
#'
#' @rdname print
#' @export
#' @method print Regions
print.Regions <- function(x, ...) {
  n_regs <- length(x@ranges)
  n_peaks <- length(unique(x@peaks))
  cat(paste0(
    "An Regions object\n", "with ", n_regs, " candidate genomic regions ",
    "in ", n_peaks, " peaks"
  ))
}

setMethod(
  "show",
  signature = "Regions",
  function(object) {
    print(object)
  }
)

.process_Network <- function(
    object,
    r_squared_threshold = 0) {
  metrics <- methods::slot(object, "metrics")
  metrics <- metrics[metrics$r_squared >= r_squared_threshold, ]
  targets <- unique(metrics$target)
  coefficients <- methods::slot(object, "coefficients")
  coefficients <- coefficients[coefficients$target %in% targets, ]

  if (is.null(coefficients) || nrow(coefficients) == 0) {
    methods::slot(object, "network") <- data.frame(
      regulator = character(),
      target = character(),
      weight = numeric(),
      mean_weight = numeric(),
      mean_corr = numeric(),
      n_regions = integer()
    )
    return(object)
  }

  coefficients_renamed <- dplyr::rename(
    coefficients,
    regulator = tf
  )
  aggregated_edges <- dplyr::group_by(
    coefficients_renamed,
    regulator,
    target
  )
  aggregated_edges <- dplyr::summarise(
    aggregated_edges,
    sum_weight = sum(coefficient),
    mean_weight = mean(coefficient),
    mean_corr = mean(corr),
    n_regions = n(),
    .groups = "drop"
  )

  aggregated_edges <- dplyr::group_by(aggregated_edges, target)
  aggregated_edges <- dplyr::mutate(
    aggregated_edges,
    weight = normalization(sum_weight, method = "unit_vector")
  )
  aggregated_edges <- dplyr::ungroup(aggregated_edges)
  aggregated_edges <- as.data.frame(aggregated_edges)

  aggregated_edges <- dplyr::select(
    aggregated_edges,
    regulator,
    target,
    weight,
    sum_weight,
    mean_weight,
    mean_corr,
    n_regions
  )

  methods::slot(object, "network") <- aggregated_edges
  return(object)
}

.process_csn <- function(
    object,
    r_squared_threshold = 0) {
  active_network <- DefaultNetwork(object)
  networks <- object@networks[[active_network]]

  if (is.null(networks) || length(networks) == 0) {
    return(object)
  }

  for (celltype in names(networks)) {
    network <- networks[[celltype]]
    if (!is.null(network) && methods::is(network, "Network")) {
      network <- .process_Network(
        network,
        r_squared_threshold = r_squared_threshold
      )
      object@networks[[active_network]][[celltype]] <- network
    }
  }

  return(object)
}

#' @param active_network Character string specifying which network to export
#' @param celltypes Character vector of cell types to export
#' @param weight_cutoff Numeric threshold for filtering edges by absolute weight
#'
#' @rdname export_csn
#' @method export_csn CSNObject
#' @export
setMethod(
  "export_csn",
  signature = "CSNObject",
  function(object,
           active_network = NULL,
           celltypes = NULL,
           weight_cutoff = NULL,
           ...) {
    active_network <- active_network %ss% DefaultNetwork(object)
    net_all <- object@networks[[active_network]]

    celltypes_all <- names(net_all)
    celltypes <- intersect(celltypes %ss% celltypes_all, celltypes_all)

    res <- lapply(
      net_all[celltypes],
      function(x) {
        export_csn(x, weight_cutoff = weight_cutoff)
      }
    ) |>
      purrr::set_names(celltypes)

    if (length(celltypes) == 1) {
      return(res[[1]])
    }

    return(res)
  }
)

#' @rdname export_csn
#' @method export_csn Network
#' @export
setMethod(
  "export_csn",
  signature = "Network",
  function(object,
           weight_cutoff = NULL,
           ...) {
    network <- methods::slot(object, "network")

    if (is.null(network)) {
      return(NULL)
    }

    network <- network_format(network, abs_weight = FALSE)

    if (!is.null(weight_cutoff)) {
      network <- network[abs(network$weight) >= weight_cutoff, ]
    }

    return(network)
  }
)
