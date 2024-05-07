#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object The input object for \code{inferCSN}
#' @param ... Arguments for other methods
#'
#' @importFrom methods as is
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict
#' @importFrom utils methods
#'
#' @return A data table of gene-gene regulatory relationship
#'
#' @rdname inferCSN
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix, verbose = TRUE)
#' head(weight_table)
#'
#' weight_table <- inferCSN(example_matrix, cores = 2)
#' head(weight_table)
#'
#' \dontrun{
#' data("promoter_regions_hg38")
#' seurat_object <- inferCSN(
#'  seurat_object,
#'  genome_info = promoter_regions_hg38
#' )
#' }
# inferCSN <- function(object, ...) {
#   UseMethod(generic = "inferCSN", object = object)
# }
setGeneric(
  "inferCSN",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "inferCSN", object = object)
  }
)

#' @title Get dimensional information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
#'
#' @rdname get.dimensional
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' dimensional_information <- get.dimensional(example_matrix)
setGeneric(
  "get.dimensional",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "get.dimensional", object = object)
  }
)

#' @title infer VECTOR
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
#'
#' @rdname inferVECTOR
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' dimensional_information <- inferVECTOR(example_matrix)
setGeneric(
  "inferVECTOR",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "inferVECTOR", object = object)
  }
)

#' @title Get dynamic genes
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
#'
#' @rdname dynamic.genes
setGeneric(
  "dynamic.genes",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "dynamic.genes", object = object)
  }
)

#' @title get_pseudotime
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
#'
#' @rdname get.pseudotime
setGeneric(
  "get.pseudotime",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "get.pseudotime", object = object)
  }
)

#' @param object The input data, a seurat object.
#' @param ... Arguments for other methods
#'
#' @rdname initiate_object
#' @export initiate_object
initiate_object <- function(object, ...) {
  UseMethod(generic = "initiate_object", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_motifs
#' @export find_motifs
find_motifs <- function(object, ...) {
  UseMethod(generic = "find_motifs", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname inferCSN2
#' @export inferCSN2
inferCSN2 <- function(object, ...) {
  UseMethod(generic = "inferCSN2", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname fit_grn_models
#' @export fit_grn_models
fit_grn_models <- function(object, ...) {
  UseMethod(generic = "fit_grn_models", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_gof
#' @export plot_gof
plot_gof <- function(object, ...) {
  UseMethod(generic = "plot_gof", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_module_metrics
#' @export plot_module_metrics
plot_module_metrics <- function(object, ...) {
  UseMethod(generic = "plot_module_metrics", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_network_graph
#' @export get_network_graph
get_network_graph <- function(object, ...) {
  UseMethod(generic = "get_network_graph", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_network_graph
#' @export plot_network_graph
plot_network_graph <- function(object, ...) {
  UseMethod(generic = "plot_network_graph", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname get_tf_network
#' @export get_tf_network
get_tf_network <- function(object, ...) {
  UseMethod(generic = "get_tf_network", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname plot_tf_network
#' @export plot_tf_network
plot_tf_network <- function(object, ...) {
  UseMethod(generic = "plot_tf_network", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkGraph
#' @export NetworkGraph
NetworkGraph <- function(object, ...) {
  UseMethod(generic = "NetworkGraph", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetGRN
#' @export GetGRN
GetGRN <- function(object, ...) {
  UseMethod(generic = "GetGRN", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetNetwork
#' @export GetNetwork
GetNetwork <- function(object, ...) {
  UseMethod(generic = "GetNetwork", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkFeatures
#' @export NetworkFeatures
NetworkFeatures <- function(object, ...) {
  UseMethod(generic = "NetworkFeatures", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkRegions
#' @export NetworkRegions
NetworkRegions <- function(object, ...) {
  UseMethod(generic = "NetworkRegions", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname Params
#' @export Params
Params <- function(object, ...) {
  UseMethod(generic = "Params", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkParams
#' @export NetworkParams
NetworkParams <- function(object, ...) {
  UseMethod(generic = "NetworkParams", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkTFs
#' @export NetworkTFs
NetworkTFs <- function(object, ...) {
  UseMethod(generic = "NetworkTFs", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname NetworkModules
#' @export NetworkModules
NetworkModules <- function(object, ...) {
  UseMethod(generic = "NetworkModules", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname DefaultNetwork
#' @export DefaultNetwork
DefaultNetwork <- function(object, ...) {
  UseMethod(generic = "DefaultNetwork", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname gof
#' @export gof
gof <- function(object, ...) {
  UseMethod(generic = "gof", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_modules
#' @export find_modules
find_modules <- function(object, ...) {
  UseMethod(generic = "find_modules", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssaySummary
#' @export GetAssaySummary
GetAssaySummary <- function(object, ...) {
  UseMethod(generic = "GetAssaySummary", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname GetAssay
#' @export GetAssay
GetAssay <- function(object, ...) {
  UseMethod(generic = "GetAssay", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname LayerData
#' @export LayerData
LayerData <- function(object, ...) {
  UseMethod(generic = "LayerData", object = object)
}

#' @param object The input data, a grn object.
#' @param ... Arguments for other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
VariableFeatures <- function(object, ...) {
  UseMethod(generic = "VariableFeatures", object = object)
}
