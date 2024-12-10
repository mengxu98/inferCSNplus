#' @title Initiate the \code{CSNObject} object
#'
#' @param object The input data, a seurat object.
#' @param ... Arguments for other methods
#'
#' @rdname initiate_object
#' @export initiate_object
setGeneric(
  "initiate_object",
  signature = "object",
  function(object, ...) {
    standardGeneric("initiate_object")
  }
)

#' @title Scan for motifs in candidate regions
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_motifs
#' @export find_motifs
setGeneric(
  "find_motifs",
  signature = "object",
  function(object, ...) {
    standardGeneric("find_motifs")
  }
)

#' @title Get network regions
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkRegions
#' @export NetworkRegions
setGeneric(
  "NetworkRegions",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkRegions")
  }
)

#' @title Get network modules
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkModules
#' @export NetworkModules
setGeneric(
  "NetworkModules",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkModules")
  }
)

#' @title Get network parameters
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkParams
#' @export NetworkParams
setGeneric(
  "NetworkParams",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkParams")
  }
)

#' @title Get network graph
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkGraph
#' @export NetworkGraph
setGeneric(
  "NetworkGraph",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkGraph")
  }
)

#' @title Get dimensional information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Embedding
#' @export
#'
#' @rdname get_embedding
setGeneric(
  "get_embedding",
  signature = "object",
  function(object, ...) {
    standardGeneric("get_embedding")
  }
)

#' @title infer VECTOR
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return vector
#' @export
#'
#' @rdname infer_vector
setGeneric(
  "infer_vector",
  signature = "object",
  function(object, ...) {
    standardGeneric("infer_vector")
  }
)

#' @title Get pseudotime information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Pseudotime information corresponding to cells
#' @export get_pseudotime
#'
#' @rdname get_pseudotime
setGeneric(
  "get_pseudotime",
  signature = "object",
  function(object, ...) {
    standardGeneric("get_pseudotime")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_gof
#' @export plot_gof
setGeneric(
  "plot_gof",
  signature = "object",
  function(object, ...) {
    standardGeneric("plot_gof")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_module_metrics
#' @export plot_module_metrics
setGeneric(
  "plot_module_metrics",
  signature = "object",
  function(object, ...) {
    standardGeneric("plot_module_metrics")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_network_graph
#' @export get_network_graph
setGeneric(
  "get_network_graph",
  signature = "object",
  function(object, ...) {
    standardGeneric("get_network_graph")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_network_graph
#' @export plot_network_graph
setGeneric(
  "plot_network_graph",
  signature = "object",
  function(object, ...) {
    standardGeneric("plot_network_graph")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_tf_network
#' @export get_tf_network
setGeneric(
  "get_tf_network",
  signature = "object",
  function(object, ...) {
    standardGeneric("get_tf_network")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_tf_network
#' @export plot_tf_network
setGeneric(
  "plot_tf_network",
  signature = "object",
  function(object, ...) {
    standardGeneric("plot_tf_network")
  }
)

#' @title Get network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetNetwork
#' @export GetNetwork
setGeneric(
  "GetNetwork",
  signature = "object",
  function(object, ...) {
    standardGeneric("GetNetwork")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname Params
#' @export Params
setGeneric(
  "Params",
  signature = "object",
  function(object, ...) {
    standardGeneric("Params")
  }
)

#' @title Get network TFs
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkTFs
#' @export NetworkTFs
setGeneric(
  "NetworkTFs",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkTFs")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname DefaultNetwork
#' @export DefaultNetwork
setGeneric(
  "DefaultNetwork",
  signature = "object",
  function(object, ...) {
    standardGeneric("DefaultNetwork")
  }
)
#' @title Get metrics
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname metrics
#' @export metrics
setGeneric(
  "metrics",
  signature = "object",
  function(object, ...) {
    standardGeneric("metrics")
  }
)

#' @title Find TF modules in regulatory network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_modules
#' @export find_modules
setGeneric(
  "find_modules",
  signature = "object",
  function(object, ...) {
    standardGeneric("find_modules")
  }
)

#' @title Get summary of seurat assay
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssaySummary
#' @export GetAssaySummary
setGeneric(
  "GetAssaySummary",
  signature = "object",
  function(object, ...) {
    standardGeneric("GetAssaySummary")
  }
)

#' @title Get seurat assay
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssay
#' @export GetAssay
setGeneric(
  "GetAssay",
  signature = "object",
  function(object, ...) {
    standardGeneric("GetAssay")
  }
)

#' @title Get layer data from CSNObject
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname LayerData
#' @export LayerData
setGeneric(
  "LayerData",
  signature = "object",
  function(object, ...) {
    standardGeneric("LayerData")
  }
)

#' @title Get variable features from CSNObject
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
setGeneric(
  "VariableFeatures",
  signature = "object",
  function(object, ...) {
    standardGeneric("VariableFeatures")
  }
)

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
