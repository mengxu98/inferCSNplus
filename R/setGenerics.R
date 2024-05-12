#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object The input object for \code{inferCSN}
#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0" and "L0L2".
#' For high-dimensional and sparse data, such as single-cell sequencing data, "L0L2" is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param cross_validation Check whether cross validation is used.
#' @param n_folds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param k_folds The number of folds for sample split.
#' @param r_threshold r_threshold.
#' @param regulators Regulator genes.
#' @param targets Target genes.
#' @param regulators_num The number of non-zore coef, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param verbose Logical. Display messages. Set verbose to '2' to print errors for all model fits.
#' @param cores CPU cores. Setting to parallelize the computation with \code{\link[foreach]{foreach}}.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A data table of gene-gene regulatory relationship
#'
#' @rdname inferCSN
#' @export
setGeneric(
  name = "inferCSN",
  signature = c("object"),
  def = function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    cores = 1,
    verbose = FALSE,
    ...) {
    UseMethod(
      generic = "inferCSN",
      object = object
    )
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

#' @title Get pseudotime information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export get.pseudotime
#'
#' @rdname get.pseudotime
setGeneric(
  "get.pseudotime",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "get.pseudotime", object = object)
  }
)

#' @title Initiate the \code{RegulatoryNetwork} object
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

#' @title Fit models for gene expression
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname fit_models
#' @export fit_models
fit_models <- function(object, ...) {
  UseMethod(generic = "fit_models", object = object)
}

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
#' @rdname GetGRN
#' @export GetGRN
GetGRN <- function(object, ...) {
  UseMethod(generic = "GetGRN", object = object)
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
#' @rdname NetworkFeatures
#' @export NetworkFeatures
NetworkFeatures <- function(object, ...) {
  UseMethod(generic = "NetworkFeatures", object = object)
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
#' @rdname gof
#' @export gof
gof <- function(object, ...) {
  UseMethod(generic = "gof", object = object)
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
