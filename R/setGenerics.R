#' Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object The object for inferCSN
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
#' weightDT <- inferCSN(example_matrix, verbose = TRUE)
#' head(weightDT)
#'
#' weightDT <- inferCSN(example_matrix, verbose = TRUE, cores = 2)
#' head(weightDT)
#'
#' \dontrun{
#' data("promoter_regions_hg38")
#' weightDT <- inferCSN(seurat_object, enome_info <- promoter_regions_hg38)
#' }
setGeneric("inferCSN",
           signature = "object",
           function(object, ...) {
             UseMethod(generic = "inferCSN", object = object)
           })
