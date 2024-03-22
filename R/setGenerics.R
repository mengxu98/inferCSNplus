#' Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object The input object for \code{inferCSN}
#' @param ... Arguments for other methods
#'
#' @import Matrix
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
#' weight_table <- inferCSN(example_matrix, verbose = TRUE, cores = 2)
#' head(weight_table)
#'
#' \dontrun{
#' data("promoter_regions_hg38")
#' seurat_object <- inferCSN(seurat_object, genome_info = promoter_regions_hg38)
#' }
# inferCSN <-  function(object, ...) {
#   UseMethod(generic = "inferCSN", object = object)
# }
setGeneric(
  "inferCSN",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "inferCSN", object = object)
  }
)

#' Get dimensional information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
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

#' Get dimensional information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
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

#' Get dimensional information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Dimensional information
#' @export
setGeneric(
  "dynamic.genes",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "dynamic.genes", object = object)
  }
)
