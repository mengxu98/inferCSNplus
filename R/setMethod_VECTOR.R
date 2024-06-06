#' @param object The information of PCA and UMAP.
#' @param plot Logical value.
#'
#' @return vector result
#' @export
#'
#' @method inferVECTOR default
#'
#' @rdname inferVECTOR
inferVECTOR.default <- function(
    object,
    plot = TRUE,
    ...) {
  return(inferVECTOR(object))
}

#' @method inferVECTOR VECTOR_dimension
#'
#' @rdname inferVECTOR
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' dimensional_information <- get.dimensional(example_matrix)
#' vector_result <- inferVECTOR(dimensional_information)
inferVECTOR.VECTOR_dimension <- function(
    object,
    plot = TRUE,
    ...) {
  PCA <- object[["PCA"]]
  VEC <- object[["UMAP"]]

  graphics::par(mfrow = c(2, 3), mai = c(0.2, 0.2, 0.1, 0.1))

  PCA <- vector.rankPCA(PCA)

  # Define pixel
  OUT <- vector.buildGrid(VEC, N = 30, plot = plot)

  # Build network
  OUT <- vector.buildNet(OUT, CUT = 1, plot = plot)

  # Calculate Quantile Polarization (QP) score
  OUT <- vector.getValue(OUT, PCA, plot = plot)

  # Get pixel's QP score
  OUT <- vector.gridValue(OUT, plot = plot)

  # Find starting point
  OUT <- vector.autoCenter(OUT, UP = 0.9, plot = plot)

  # Infer vector
  OUT <- vector.drawArrow(
    OUT,
    P = 0.9,
    plot = plot,
    COL = OUT$COL,
    plot.SUMMIT = plot
  )

  graphics::par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

  graphics::mtext("UMAP1", side = 1)
  graphics::mtext("UMAP2", side = 2)

  return(OUT)
}

#' @param filter Logical value, default set to \code{TRUE},
#' will filter cells though pseudotime.
#'
#' @method inferVECTOR matrix
#'
#' @rdname inferVECTOR
#'
#' @examples
#' data("example_matrix")
#' vector_result <- inferVECTOR(
#'   example_matrix,
#'   plot = FALSE
#' )
inferVECTOR.matrix <- function(
    object,
    filter = TRUE,
    plot = TRUE,
    ...) {
  dimension <- get.dimensional(object)
  result <- inferVECTOR(dimension)
  selected <- which(result$P.PS != "NA")
  object <- object[selected, ]

  if (is.null(rownames(object))) {
    rownames(object) <- paste0(rep("cell_", nrow(object)), 1:nrow(object))
  }
  pseudotime <- result$P.PS[selected]
  pseudotime <- normalization(pseudotime)

  object <- list(
    matrix = object,
    pseudotime = data.frame(
      cells = rownames(object),
      pseudotime = pseudotime
    )
  )

  return(object)
}

#' @method inferVECTOR Seurat
#'
#' @rdname inferVECTOR
inferVECTOR.Seurat <- function(
    object,
    filter = TRUE,
    plot = TRUE,
    ...) {
  dimension <- get.dimensional(
    object
  )
  result <- inferVECTOR(dimension)
  selected <- which(result$P.PS != "NA")
  pseudotime <- result$P.PS[selected]
  pseudotime <- normalization(pseudotime)
  object <- object[, selected]
  object@meta.data$pseudotime <- pseudotime

  return(object)
}
