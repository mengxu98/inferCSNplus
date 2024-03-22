#' The Modules class
#'
#' The Modules object stores the TF modules extracted from the inferred network.
#'
#' @slot matrix Filtered matrix.
#' @slot pseduotime Pseduotime of cells
#'
#' @name VECTOR-class
#' @rdname VECTOR-class
#' @exportClass VECTOR
setClass(
  Class = "VECTOR",
  slots = list(
    matrix = "matrix",
    pseudotime = "data.frame"
  )
)
