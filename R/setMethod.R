#' @include setClass.R

#' @title Get network
#'
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname GetNetwork
#' @method GetNetwork CSNObject
#' @export
GetNetwork.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    celltypes = NULL,
    ...) {
  if (!is.null(celltypes)) {
    return(lapply(celltypes, function(ct) object@networks[[network]][[ct]]))
  }
  return(object@networks[[network]])
}

#' @title Get network TFs
#'
#' @rdname NetworkTFs
#' @method NetworkTFs CSNObject
#' @export
NetworkTFs.CSNObject <- function(object, ...) {
  return(object@regions@motifs2tfs)
}

#' @title Get network regions
#'
#' @rdname NetworkRegions
#' @method NetworkRegions CSNObject
#' @export
NetworkRegions.CSNObject <- function(object, ...) {
  return(object@regions)
}

#' @title Get TF modules
#'
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkModules
#' @method NetworkModules CSNObject
#' @export
NetworkModules.CSNObject <- function(
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
    return(lapply(networks, function(net) net@modules))
  }
  return(lapply(networks, function(net) net@modules))
}

#' @rdname NetworkModules
#' @method NetworkModules Network
#' @export
NetworkModules.Network <- function(object, ...) {
  return(object@modules)
}

#' @title Get network parameters
#'
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkParams
#' @method NetworkParams CSNObject
#' @export
NetworkParams.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    celltypes = NULL,
    ...) {
  networks <- GetNetwork(object, network = network, celltypes = celltypes)
  if (is.null(celltypes)) {
    return(lapply(networks, function(net) net@params))
  }
  return(lapply(networks, function(net) net@params))
}

#' @rdname NetworkParams
#' @method NetworkParams Network
#' @export
NetworkParams.Network <- function(object, ...) {
  return(object@params)
}

#' @title Get network parameters
#'
#' @param network network
#' @param graph graph
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname NetworkGraph
#' @method NetworkGraph CSNObject
#' @export
NetworkGraph.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    graph = "module_graph",
    celltypes = NULL,
    ...) {
  networks <- GetNetwork(object, network = network, celltypes = celltypes)
  if (is.null(celltypes)) {
    return(lapply(networks, function(net) NetworkGraph(net, graph = graph)))
  }
  return(lapply(networks, function(net) NetworkGraph(net, graph = graph)))
}

#' @param graph graph
#'
#' @rdname NetworkGraph
#' @method NetworkGraph Network
#' @export
NetworkGraph.Network <- function(
    object,
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

#' @title Get active network
#'
#' @rdname DefaultNetwork
#' @method DefaultNetwork CSNObject
#' @export
DefaultNetwork.CSNObject <- function(object, ...) {
  return(object@active_network)
}

#' @title Get fitted coefficients
#'
#' @param object The csn object
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#' @param ... other parameters
#'
#' @rdname coef
#' @method coef CSNObject
#' @export
coef.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    celltypes = NULL,
    ...) {
  networks <- GetNetwork(object, network = network, celltypes = celltypes)
  if (is.null(celltypes)) {
    return(lapply(networks, function(net) net@coefficients))
  }
  return(lapply(networks, function(net) net@coefficients))
}

#' @rdname coef
#' @method coef Network
#' @export
coef.Network <- function(object, ...) {
  return(object@coefficients)
}

#' @title Get goodness-of-fit info
#'
#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#'
#' @rdname metrics
#' @method metrics CSNObject
#' @export
metrics.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    celltypes = NULL,
    ...) {
  networks <- GetNetwork(object, network = network, celltypes = celltypes)
  if (is.null(celltypes)) {
    return(lapply(networks, function(net) net@fit))
  }
  return(lapply(networks, function(net) net@fit))
}

#' @rdname metrics
#' @method metrics Network
#' @export
metrics.Network <- function(
    object,
    celltypes = NULL,
    ...) {
  return(object@fit)
}

#' @title Get GRN inference parameters
#'
#' @rdname Params
#' @method Params CSNObject
#' @export
Params.CSNObject <- function(object, ...) {
  return(object@params)
}

#' @title Get summary of seurat assay
#'
#' @param group_name group_name
#' @param assay assay
#' @param verbose verbose
#'
#' @rdname GetAssaySummary
#' @method GetAssaySummary Seurat
#' @export
GetAssaySummary.Seurat <- function(
    object,
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

#' @title Get summary of seurat assay from CSNObject
#'
#' @param group_name group_name
#' @param assay assay
#' @param verbose verbose
#'
#' @rdname GetAssaySummary
#' @method GetAssaySummary CSNObject
#' @export
GetAssaySummary.CSNObject <- function(
    object,
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

#' @title Get Seurat assay from CSNObject
#'
#' @param assay assay
#'
#' @rdname GetAssay
#' @method GetAssay CSNObject
#' @export
GetAssay.CSNObject <- function(
    object,
    assay = NULL,
    ...) {
  return(
    Seurat::GetAssay(
      object@data,
      assay = assay
    )
  )
}

#' @title Get layer data from CSNObject
#'
#' @param ... other parameters
#'
#' @rdname LayerData
#' @method LayerData CSNObject
#' @export
LayerData.CSNObject <- function(object, ...) {
  return(SeuratObject::LayerData(object@data, ...))
}

#' @title Get variable features from CSNObject
#'
#' @rdname VariableFeatures
#' @method VariableFeatures CSNObject
#' @export
VariableFeatures.CSNObject <- function(object, ...) {
  return(Seurat::VariableFeatures(object@data, ...))
}

#' @title Print Network objects
#'
#' @param x x
#' @param ... other parameters
#'
#' @rdname print
#' @export
#' @method print Network
print.Network <- function(x, ...) {
  if (nrow(NetworkModules(x)@meta) == 0) {
    n_genes <- length(unique(coef(x)$target))
    n_tfs <- length(unique(coef(x)$tf))
  } else {
    n_genes <- length(unique(NetworkModules(x)@meta$target))
    n_tfs <- length(unique(NetworkModules(x)@meta$tf))
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

# setMethod(
#   "show",
#   signature = "Regions",
#   function(object) {
#     print(object)
#   }
# )


#' @title Get network data
#' @rdname summary_csn
#' @method summary_csn CSNObject
#' @export
summary_csn.CSNObject <- function(object, ...) {
  networks <- object@networks

  tfs <- unique(unlist(lapply(networks, function(network) {
    if (!is.null(network) && methods::is(network, "Network")) {
      network_df <- methods::slot(network, "network")
      if (!is.null(network_df) && nrow(network_df) > 0) {
        unique(network_df$regulator)
      }
    }
  })))

  network_stats <- lapply(names(networks), function(celltype) {
    network <- networks[[celltype]]
    if (!is.null(network) && methods::is(network, "Network")) {
      network_df <- methods::slot(network, "network")
      if (!is.null(network_df) && nrow(network_df) > 0) {
        data.frame(
          celltype = celltype,
          edges = nrow(network_df),
          targets = length(unique(network_df$target)),
          regulators = length(unique(network_df$regulator))
        )
      }
    }
  })
  network_stats <- do.call(rbind, network_stats)

  msg <- sprintf("A CSNObject object based on %d transcription factors\n\n", length(tfs))
  msg <- paste0(msg, sprintf(
    "%d inferred networks: %s\n\n",
    length(names(networks)),
    paste(names(networks), collapse = ", ")
  ))

  if (!is.null(network_stats)) {
    msg <- paste0(msg, "Network statistics by celltype:\n")
    for (i in 1:nrow(network_stats)) {
      msg <- paste0(msg, sprintf(
        "%s: %d edges (%d regulators -> %d targets)\n",
        network_stats$celltype[i],
        network_stats$edges[i],
        network_stats$regulators[i],
        network_stats$targets[i]
      ))
    }
  }

  result <- list(
    networks = networks,
    tfs = tfs,
    network_stats = network_stats
  )

  class(result) <- "CSNSummary"

  return(result)
}

#' @export
print.CSNSummary <- function(x, ...) {
  cat(sprintf("$networks\n%s\n\n", x$msg))
  cat("$tfs\n")
  print(x$tfs)
  cat("\n$network_stats\n")
  print(x$network_stats)
}

#' @title Get regulatory genes
#'
#' @param object A CSNObject object
#'
#' @rdname get_attribute
#' @export
setGeneric(
  "get_attribute",
  function(object, ...) {
    standardGeneric("get_attribute")
  }
)

#' @md
#' @param celltypes A character vector specifying the celltypes to get attributes for.
#' If \code{NULL}, all celltypes are returned.
#' @param active_network A character string specifying the active network to get attributes for.
#' @param attribute A character string specifying the attribute to get.
#' This can take any of the following choices:
#'
#' * *`genes`* - Gene names
#'
#' * *`tfs`* - Transcription factors
#'
#' * *`peaks`* - Peaks
#'
#' * *`regulators`* - Regulators
#'
#' * *`targets`* - Targets
#'
#' * *`cells`* - Cells
#'
#' * *`modules`* - Modules
#'
#' * *`coefficients`* - Coefficients
#'
#' @param ... Additional arguments
#'
#' @return A character vector of regulatory genes
#'
#' @rdname get_attribute
#' @method get_attribute CSNObject
#'
#' @export
setMethod(
  "get_attribute",
  "CSNObject",
  function(object,
           celltypes = NULL,
           active_network = NULL,
           attribute = c(
             "genes",
             "tfs",
             "peaks",
             "regulators",
             "targets",
             "celltypes",
             "cells",
             "modules",
             "coefficients"
           ),
           ...) {
    if (is.null(celltypes)) {
      celltypes <- object@metadata$celltypes
    }
    if (is.null(active_network)) {
      active_network <- DefaultNetwork(object)
    }
    attribute <- match.arg(attribute)
    attributes <- switch(
      EXPR = attribute,
      "genes" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$genes$gene
        }
      ),
      "peaks" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$peaks$peak
        }
      ),
      "regulators" = lapply(
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]$regulators |>
            unique()
        }
      ),
      "targets" = lapply(
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]$targets |>
            unique()
        }
      ),
      "cells" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$cells
        }
      ),
      "modules" = lapply( # TODO: check if this is correct
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]$modules
        }
      ),
      "coefficients" = lapply( # TODO: check if this is correct
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]$coefficients
        }
      )
    )
    if (attribute == "tfs") {
      attributes <- object@metadata$tfs
      return(attributes)
    }

    attributes <- purrr::set_names(attributes, celltypes)
    if (length(attributes) == 1) {
      return(attributes[[1]])
    }

    return(attributes)
  }
)

.process_csn <- function(object) {
  active_network <- DefaultNetwork(object)
  networks <- object@networks[[active_network]]

  if (is.null(networks) || length(networks) == 0) {
    return(object)
  }

  for (celltype in names(networks)) {
    network <- networks[[celltype]]
    if (!is.null(network) && methods::is(network, "Network")) {
      coefficients <- methods::slot(network, "coefficients")

      if (is.null(coefficients) || nrow(coefficients) == 0) {
        methods::slot(network, "network") <- data.frame(
          regulator = character(),
          target = character(),
          weight = numeric(),
          mean_weight = numeric(),
          mean_corr = numeric(),
          n_regions = integer()
        )
        object@networks[[active_network]][[celltype]] <- network
        next
      }

      # First calculate aggregated edges with sum_weight
      aggregated_edges <- coefficients |>
        dplyr::rename(regulator = tf) |>
        dplyr::group_by(regulator, target) |>
        dplyr::summarise(
          sum_weight = sum(coefficient),
          mean_weight = mean(coefficient),
          mean_corr = mean(corr),
          n_regions = n()
        ) |>
        dplyr::ungroup()

      # Then normalize sum_weight for each target and save as weight
      aggregated_edges <- aggregated_edges |>
        dplyr::group_by(target) |>
        dplyr::mutate(
          weight = normalization(
            sum_weight,
            method = "unit_vector"
          )
        ) |>
        dplyr::ungroup() |>
        as.data.frame()

      aggregated_edges <- aggregated_edges |>
        dplyr::select(
          regulator,
          target,
          weight,
          sum_weight,
          mean_weight,
          mean_corr,
          n_regions
        )

      methods::slot(network, "network") <- aggregated_edges
      object@networks[[active_network]][[celltype]] <- network
    }
  }

  return(object)
}

#' @title Export network from CSN object
#' @description Export network data from a CSN object with optional filtering
#'
#' @param object A CSNObject object
#' @param ... Additional arguments
#'
#' @return A list of network data frames by cell type
#'
#' @export
setGeneric(
  "export_csn",
  function(object, ...) {
    standardGeneric("export_csn")
  }
)

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
    active_network <- active_network %s% DefaultNetwork(object)
    net_all <- object@networks[[active_network]]

    celltypes_all <- names(net_all)
    celltypes <- intersect(celltypes %s% celltypes_all, celltypes_all)

    res <- lapply(
      net_all[celltypes],
      function(x) {
        network <- methods::slot(x, "network") |>
          network_format(abs_weight = FALSE)

        if (is.null(network)) {
          return(NULL)
        }

        if (!is.null(weight_cutoff)) {
          network <- network[abs(network$weight) >= weight_cutoff, ]
        }

        return(network)
      }
    ) |>
      purrr::set_names(celltypes)

    if (length(celltypes) == 1) {
      return(res[[1]])
    }

    return(res)
  }
)
