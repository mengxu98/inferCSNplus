#' @title Initiate the \code{CSNObject} object
#'
#' @param object The input data, a seurat object.
#' @param ... Arguments for other methods
#'
#' @rdname initiate_object
#' @export initiate_object
initiate_object <- function(object, ...) {
  UseMethod(generic = "initiate_object", object = object)
}

#' @title Scan for motifs in candidate regions
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_motifs
#' @export find_motifs
find_motifs <- function(object, ...) {
  UseMethod(generic = "find_motifs", object = object)
}

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
    UseMethod(generic = "get_embedding", object = object)
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
    UseMethod(generic = "infer_vector", object = object)
  }
)

#' @title Get dynamic genes
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return A new object with dynamic genes
#' @export
#'
#' @rdname dynamic_genes
setGeneric(
  "dynamic_genes",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "dynamic_genes", object = object)
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
    UseMethod(generic = "get_pseudotime", object = object)
  }
)



#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_gof
#' @export plot_gof
plot_gof <- function(object, ...) {
  UseMethod(generic = "plot_gof", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_module_metrics
#' @export plot_module_metrics
plot_module_metrics <- function(object, ...) {
  UseMethod(generic = "plot_module_metrics", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_network_graph
#' @export get_network_graph
get_network_graph <- function(object, ...) {
  UseMethod(generic = "get_network_graph", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_network_graph
#' @export plot_network_graph
plot_network_graph <- function(object, ...) {
  UseMethod(generic = "plot_network_graph", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_tf_network
#' @export get_tf_network
get_tf_network <- function(object, ...) {
  UseMethod(generic = "get_tf_network", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_tf_network
#' @export plot_tf_network
plot_tf_network <- function(object, ...) {
  UseMethod(generic = "plot_tf_network", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkGraph
#' @export NetworkGraph
NetworkGraph <- function(object, ...) {
  UseMethod(generic = "NetworkGraph", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname summary_csn
#' @export summary_csn
summary_csn <- function(object, ...) {
  UseMethod(generic = "summary_csn", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetNetwork
#' @export GetNetwork
GetNetwork <- function(object, ...) {
  UseMethod(generic = "GetNetwork", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkRegions
#' @export NetworkRegions
NetworkRegions <- function(object, ...) {
  UseMethod(generic = "NetworkRegions", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname Params
#' @export Params
Params <- function(object, ...) {
  UseMethod(generic = "Params", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkParams
#' @export NetworkParams
NetworkParams <- function(object, ...) {
  UseMethod(generic = "NetworkParams", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkTFs
#' @export NetworkTFs
NetworkTFs <- function(object, ...) {
  UseMethod(generic = "NetworkTFs", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkModules
#' @export NetworkModules
NetworkModules <- function(object, ...) {
  UseMethod(generic = "NetworkModules", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname DefaultNetwork
#' @export DefaultNetwork
DefaultNetwork <- function(object, ...) {
  UseMethod(generic = "DefaultNetwork", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname metrics
#' @export metrics
metrics <- function(object, ...) {
  UseMethod(generic = "metrics", object = object)
}

#' @title Find TF modules in regulatory network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_modules
#' @export find_modules
find_modules <- function(object, ...) {
  UseMethod(generic = "find_modules", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssaySummary
#' @export GetAssaySummary
GetAssaySummary <- function(object, ...) {
  UseMethod(generic = "GetAssaySummary", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssay
#' @export GetAssay
GetAssay <- function(object, ...) {
  UseMethod(generic = "GetAssay", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname LayerData
#' @export LayerData
LayerData <- function(object, ...) {
  UseMethod(generic = "LayerData", object = object)
}

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
VariableFeatures <- function(object, ...) {
  UseMethod(generic = "VariableFeatures", object = object)
}
