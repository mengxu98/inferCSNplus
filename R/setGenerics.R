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
setGeneric("inferCSN",
           signature = "object",
           function(object, ...) {
             UseMethod(generic = "inferCSN", object = object)
           })

#' Get embedding information
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return Embedding information
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' embedding_information <- get.embedding(example_matrix)
setGeneric("get.embedding",
           signature = "object",
           function(object, ...) {
             UseMethod(generic = "get.embedding", object = object)
           })
