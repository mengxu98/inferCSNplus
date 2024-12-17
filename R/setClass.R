#' @title The Modules class
#' @description
#'  The Modules object stores the TF modules extracted from the inferred network.
#'
#' @slot meta A dataframe with meta data about the modules.
#' @slot features A named list with a set of fetures (genes/peaks) for each module.
#' @slot params A named list with module selection parameters.
#'
#' @name Modules-class
#' @rdname Modules-class
#' @exportClass Modules
setClass(
  Class = "Modules",
  slots = list(
    meta = "data.frame",
    features = "list",
    params = "list"
  )
)

#' @title The Network class
#' @description
#'  The Network object stores the inferred network itself,
#'  information about the fitting process as well as graph representations of the network.
#'
#' @slot data A matrix.
#' @slot modules A list TF modules.
#' @slot regulators A named list containing the transcription factors included in the network.
#' @slot targets A named list containing the target genes included in the network.
#' @slot metrics A dataframe with goodness of fit measures.
#' @slot coefficients A dataframe with the fitted coefficients.
#' @slot graphs Graphical representations of the inferred network.
#' @slot network A dataframe containing the network edges and their properties.
#' @slot params A named list with GRN inference parameters.
#'
#' @name Network-class
#' @rdname Network-class
#' @exportClass Network
setClass(
  Class = "Network",
  slots = list(
    data = "ANY",
    params = "list",
    regulators = "character",
    targets = "character",
    metrics = "data.frame",
    coefficients = "data.frame",
    network = "data.frame",
    modules = "Modules",
    graphs = "list"
  )
)

#' @title The Regions class
#' @description
#'  The Regions object stores the genomic regions that are considered by the model.
#'  It stores their genomic positions, how they map to the peaks in the Seurat object
#'  and motif matches.
#'
#' @slot motifs A \code{Motifs} object with matches of TF motifs.
#' @slot motifs2tfs tfs.
#' @slot ranges A \code{GenomicRanges} object.
#' @slot peaks A numeric vector with peak indices for each region.
#'
#' @name Regions-class
#' @rdname Regions-class
#' @exportClass Regions
setClass(
  Class = "Regions",
  slots = list(
    motifs = "ANY",
    motifs2tfs = "ANY",
    ranges = "GRanges",
    peaks = "numeric"
  )
)

#' @title The CSNObject class
#' @description
#'  The CSNObject object is an extended \code{Seurat} object
#'  for the storage and analysis of celltype-specific gene regulatory networks.
#'
#' @slot data Seurat object.
#' @slot metadata A list with metadata about the object.
#' @slot regions A \code{\linkS4class{Regions}} object containing information about
#' the genomic regions included in the network.
#' @slot networks A \code{\linkS4class{Network}} object containing the inferred regulatory
#' network and information about the model fit.
#' @slot params A list storing parameters for network inference.
#' @slot active_network A string indicating the active network.
#'
#' @name CSNObject-class
#' @rdname CSNObject-class
#' @exportClass CSNObject
setClass(
  Class = "CSNObject",
  slots = list(
    data = "Seurat",
    metadata = "list",
    regions = "Regions",
    networks = "list",
    params = "list",
    active_network = "character"
  )
)
