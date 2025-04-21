#' @title Get any attribute from a CSNObject object
#'
#' @param object A CSNObject object
#' @param ... Additional arguments
#'
#' @rdname get_attribute
#' @export
setGeneric(
  "get_attribute",
  function(object, ...) {
    standardGeneric("get_attribute")
  }
)

#' @param celltypes A character vector specifying the celltypes to get attributes for.
#' If \code{NULL}, all celltypes are returned.
#' @param active_network A character string specifying the active network to get attributes for.
#' @param attribute A character string specifying the attribute to get.
#' This can take any of the following choices:
#' \describe{
#' \item{genes}{The original gene without any filtering.}
#' \item{tfs}{The original transcription factors without any filtering.}
#' \item{peaks}{The original peaks without any filtering.}
#' \item{regulators}{The regulators of the network after inference CSN.}
#' \item{targets}{The targets of the network after inference CSN.}
#' \item{cells}{The cells of the network.}
#' \item{modules}{The modules of the network.}
#' \item{coefficients}{The coefficients of the network.}
#' }
#'
#' @return A character vector of regulatory genes
#'
#' @rdname get_attribute
#' @method get_attribute CSNObject
#'
#' @export
setMethod(
  "get_attribute",
  "CSNObject",
  function(object,
           celltypes = NULL,
           active_network = NULL,
           attribute = c(
             "genes",
             "tfs",
             "peaks",
             "regulators",
             "targets",
             "celltypes",
             "cells",
             "modules",
             "coefficients"
           ),
           ...) {
    if (is.null(celltypes)) {
      celltypes <- object@metadata$celltypes
    }
    if (is.null(active_network)) {
      active_network <- DefaultNetwork(object)
    }
    attribute <- match.arg(attribute)
    attributes <- switch(
      EXPR = attribute,
      "genes" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$genes$gene
        }
      ),
      "peaks" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$peaks$peak
        }
      ),
      "regulators" = {
        networks <- GetNetwork(
          object,
          network = active_network,
          celltypes = celltypes
        )
        lapply(
          networks,
          function(net) {
            coeffs <- methods::slot(net, "coefficients")
            modules_obj <- methods::slot(net, "modules")
            modules_meta <- methods::slot(modules_obj, "meta")
            if (nrow(modules_meta) == 0) {
              unique(coeffs$tf)
            } else {
              unique(modules_meta$tf)
            }
          }
        )
      },
      "targets" = {
        networks <- GetNetwork(
          object,
          network = active_network,
          celltypes = celltypes
        )
        lapply(
          networks,
          function(net) {
            coeffs <- methods::slot(net, "coefficients")
            modules_obj <- methods::slot(net, "modules")
            modules_meta <- methods::slot(modules_obj, "meta")
            if (nrow(modules_meta) == 0) {
              unique(coeffs$target)
            } else {
              unique(modules_meta$target)
            }
          }
        )
      },
      "cells" = lapply(
        celltypes,
        function(c) {
          object@metadata$attributes[[c]]$cells
        }
      ),
      "modules" = lapply( # TODO: check if this is correct
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]@modules
        }
      ),
      "coefficients" = lapply( # TODO: check if this is correct
        celltypes,
        function(c) {
          object@networks[[active_network]][[c]]@coefficients
        }
      )
    )
    if (attribute == "tfs") {
      attributes <- object@metadata$tfs
      return(attributes)
    }
    if (attribute == "celltypes") {
      attributes <- object@metadata$celltypes
      return(attributes)
    }

    attributes <- purrr::set_names(attributes, celltypes)
    if (length(attributes) == 1) {
      return(attributes[[1]])
    }

    return(attributes)
  }
)
