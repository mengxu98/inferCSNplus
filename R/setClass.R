#' @title The Modules class
#' @description
#'  The Modules object stores the TF modules extracted from the inferred network..
#'
#' @slot meta A dataframe with meta data about the modules.
#' @slot features A named list with a set of fetures (genes/peaks) for each module.
#' @slot params A named list with module selection parameters.
#'
#' @name Modules-class
#' @rdname Modules-class
#' @exportClass Modules
Modules <- setClass(
  Class = "Modules",
  slots = list(
    meta = "data.frame",
    features = "list",
    params = "list"
  )
)

#' @title The Network class
#'@description
#'  The Network object stores the inferred network itself, information about the fitting
#' process as well as graph representations of the network.
#'
#' @slot modules A list TF modules.
#' @slot features A named list containing the transcription factors and
#' target genes included in the network.
#' @slot fit A dataframe with goodness of fit measures.
#' @slot coefs A dataframe with the fitted coefficients.
#' @slot graphs Graphical representations of the inferred network.
#' @slot params A named list with GRN inference parameters.
#'
#' @name Network-class
#' @rdname Network-class
#' @exportClass Network
Network <- setClass(
  Class = "Network",
  slots = list(
    modules = "Modules",
    features = "character",
    fit = "data.frame",
    coefs = "data.frame",
    graphs = "list",
    params = "list"
  )
)

#' @title The Regions class
#' @description
#'  The Regions object stores the genomic regions that are considered by the model.
#' It stores their genomic positions, how they map to the peaks in the Seurat object
#' and motif matches.
#'
#' @slot motifs A \code{Motifs} object with matches of TF motifs.
#' @slot tfs tfs.
#' @slot ranges A \code{GenomicRanges} object.
#' @slot peaks A numeric vector with peak indices for each region.
#'
#' @name Regions-class
#' @rdname Regions-class
#' @exportClass Regions
Regions <- setClass(
  Class = "Regions",
  slots = list(
    motifs = "ANY",
    tfs = "ANY",
    ranges = "GRanges",
    peaks = "numeric"
  )
)


#' @title The RegularotyNetwork class
#' @description
#'  The RegularotyNetwork object is the core data structure in Pando.
#' It stores all data necessary for network inference and analysis
#' that is not provided by Seurat.
#'
#' @slot regions A \code{\link{Regions}} object containing information about
#' the genomic regions included in the network.
#' @slot networks A \code{\link{Network}} object containing the inferred regulatory
#' network and information about the model fit.
#' @slot params A list storing parameters for GRN inference.
#' @slot active_network A string indicating the active network.
#'
#' @name RegulatoryNetwork-class
#' @rdname RegulatoryNetwork-class
#' @exportClass RegulatoryNetwork
RegulatoryNetwork <- setClass(
  Class = "RegulatoryNetwork",
  slots = list(
    regions = "Regions",
    networks = "list",
    params = "list",
    active_network = "character"
  )
)


#' @title The CSNObject class
#' @description
#'  The CSNObject object is an extended \code{Seurat} object
#' for the storage and analysis of Regulatory network data.
#'
#' @slot csn A named list containing \code{RegulatoryNetwork} objects with inferred networks.
#' @slot data Seurat object.
#'
#' @name CSNObject-class
#' @rdname CSNObject-class
#' @exportClass CSNObject
CSNObject <- setClass(
  Class = "CSNObject",
  slots = list(
    "csn" = "RegulatoryNetwork",
    "data" = "Seurat"
  )
)

#' @title The CSNObjectList class
#' @description
#'  The CSNObjectList object is an extended \code{Seurat} object
#' for the storage and analysis of Regulatory network data.
#'
#' @slot data Seurat object list.
#'
#' @name CSNObjectList-class
#' @rdname CSNObjectList-class
#' @exportClass CSNObjectList
CSNObjectList <- setClass(
  Class = "CSNObjectList",
  slots = list(
    "data" = "list"
  )
)

#' @title Get network
#'
#' @param network network
#'
#' @rdname GetNetwork
#' @method GetNetwork CSNObject
#' @export
GetNetwork.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object@csn, network = network))
}

#' @param network network
#'
#' @rdname GetNetwork
#' @method GetNetwork RegulatoryNetwork
#' @export
GetNetwork.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  if (length(network) == 0) {
    stop(paste0("Network not found, try running `infer_network()` first."))
  }
  if (!network %in% names(object@networks)) {
    stop(paste0('The requested network "', network, '" does not exist.'))
  }
  return(object@networks[[network]])
}

#' @title Get network features
#'
#' @param network network
#'
#' @rdname NetworkFeatures
#' @method NetworkFeatures CSNObject
#' @export
NetworkFeatures.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@features)
}

#' @rdname NetworkFeatures
#' @method NetworkFeatures RegulatoryNetwork
#' @export
NetworkFeatures.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@features)
}


#' @title Get network TFs
#' 
#' @rdname NetworkTFs
#' @method NetworkTFs CSNObject
#' @export
NetworkTFs.CSNObject <- function(object, ...) {
  return(object@csn@regions@tfs)
}


#' @rdname NetworkTFs
#' @method NetworkTFs RegulatoryNetwork
#' @export
NetworkTFs.RegulatoryNetwork <- function(object, ...) {
  return(object@regions@tfs)
}


#' @title Get network regions
#'
#' @rdname NetworkRegions
#' @method NetworkRegions CSNObject
#' @export
NetworkRegions.CSNObject <- function(object, ...) {
  return(object@csn@regions)
}


#' @rdname NetworkRegions
#' @method NetworkRegions RegulatoryNetwork
#' @export
NetworkRegions.RegulatoryNetwork <- function(object, ...) {
  return(object@regions)
}


#' @title Get network data
#' @rdname GetGRN
#' @method GetGRN CSNObject
#' @export
GetGRN.CSNObject <- function(object, ...) {
  return(object@csn)
}

#' @title Get TF modules
#'
#' @param network network
#'
#' @rdname NetworkModules
#' @method NetworkModules CSNObject
#' @export
NetworkModules.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@modules)
}

#' @param network network
#'
#' @rdname NetworkModules
#' @method NetworkModules RegulatoryNetwork
#' @export
NetworkModules.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@modules)
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
#'
#' @rdname NetworkParams
#' @method NetworkParams CSNObject
#' @export
NetworkParams.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@params)
}

#' @param network network
#'
#' @rdname NetworkParams
#' @method NetworkParams RegulatoryNetwork
#' @export
NetworkParams.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@params)
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
#'
#' @rdname NetworkGraph
#' @method NetworkGraph CSNObject
#' @export
NetworkGraph.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    graph = "module_graph",
    ...) {
  return(
    NetworkGraph(
      GetNetwork(object, network = network),
      graph = graph
    )
  )
}

#' @param network network
#' @param graph graph
#'
#' @rdname NetworkGraph
#'
#' @method NetworkGraph RegulatoryNetwork
#' @export
NetworkGraph.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    graph = "module_graph",
    ...) {
  return(
    NetworkGraph(
      GetNetwork(object, network = network),
      graph = graph
    )
  )
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
      paste0("The requested graph '", graph, "' does not exist. \nTry (re-)running `get_network_graph()`.")
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
  return(DefaultNetwork(GetGRN(object)))
}

#' @rdname DefaultNetwork
#' @method DefaultNetwork RegulatoryNetwork
#' @export
DefaultNetwork.RegulatoryNetwork <- function(object, ...) {
  return(object@active_network)
}

#' @title Get fitted coefficients
#'
#' @param object The csn object
#' @param network network
#' @param ... other parameters
#'
#' @rdname coef
#' @method coef CSNObject
#' @export
coef.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@coefs)
}

#' @param network network
#'
#' @rdname coef
#' @method coef RegulatoryNetwork
#' @export
coef.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@coefs)
}

#' @rdname coef
#' @method coef Network
#' @export
coef.Network <- function(object, ...) {
  return(object@coefs)
}

#' @title Get goodness-of-fit info
#'
#' @param network network
#'
#' @rdname gof
#' @method gof CSNObject
#' @export
gof.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@fit)
}

#' @rdname gof
#' @method gof RegulatoryNetwork
#' @export
gof.RegulatoryNetwork <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(GetNetwork(object, network = network)@fit)
}

#' @rdname gof
#' @method gof Network
#' @export
gof.Network <- function(
    object,
    network = DefaultNetwork(object),
    ...) {
  return(object@fit)
}

#' @title Get GRN inference parameters
#'
#' @rdname Params
#' @method Params CSNObject
#' @export
Params.CSNObject <- function(object, ...) {
  return(object@csn@params)
}

#' @rdname Params
#' @method Params RegulatoryNetwork
#' @export
Params.RegulatoryNetwork <- function(object, ...) {
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
  smry <- Seurat::Misc(object[[assay]])$summary[[group_name]]
  if (is.null(smry)) {
    log_message('Summary of "', group_name, '" does not yet exist.', verbose = verbose)
    log_message("Summarizing.", verbose = verbose)
    object <- aggregate_assay(object, assay = assay, group_name = group_name)
    smry <- GetAssaySummary(object, assay = assay, group_name = group_name, verbose = verbose)
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
  return(GetAssaySummary(object@data, group_name, assay = NULL, verbose = TRUE))
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
  return(Seurat::GetAssay(object@data, assay = assay))
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

#' @title Print RegulatoryNetwork objects
#'
#' @param x A csn object
#' @param ... other parameters
#'
#' @rdname print
#' @export
#' @method print RegulatoryNetwork
print.RegulatoryNetwork <- function(x, ...) {
  n_tfs <- ncol(NetworkTFs(x))
  if (is.null(n_tfs)) {
    tf_string <- "\nCandidate regions have not been scanned for motifs"
  } else {
    tf_string <- paste0("based on ", n_tfs, " transcription factors\n")
  }
  n_nets <- length(x@networks)
  net_names <- names(x@networks)
  if (n_nets == 0) {
    conn_string <- "\nNo network has been inferred\n"
  } else if (n_nets == 1) {
    conn_string <- paste0(n_nets, " inferred network: ", net_names, "\n")
  } else {
    conn_string <- paste0(
      n_nets, " inferred networks: ",
      paste(net_names, collapse = ", "), "\n"
    )
  }
  cat(paste0(
    "A RegulatoryNetwork object ", tf_string, "\n",
    conn_string
  ))
}

setMethod(
  "show",
  signature = "RegulatoryNetwork",
  function(object) {
    print(object)
  }
)

#' @title Print Network objects
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

setMethod(
  "show",
  signature = "Regions",
  function(object) {
    print(object)
  }
)
