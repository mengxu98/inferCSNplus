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
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE)
#' head(weightDT)
#'
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE, cores = 2)
#' head(weightDT)
#'
inferCSN <- function(object, ...) {
  UseMethod(generic = 'inferCSN', object = object)
}
