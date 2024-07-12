#' @param object The information of pca and UMAP.
#' @param plot Logical value.
#'
#' @export
#'
#' @method infer_vector default
#'
#' @rdname infer_vector
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' dimensional_information <- get_embedding(example_matrix)
#' vector_result <- infer_vector(dimensional_information)
#' }
setMethod(
  f = "infer_vector",
  signature = signature(object = "list"),
  definition = function(
    object,
    plot = TRUE,
    ...) {
    pca <- object[["pca"]]
    umap <- object[["umap"]]

    graphics::par(mfrow = c(2, 3), mai = c(0.2, 0.2, 0.1, 0.1))

    pca <- vector_rank_pca(pca)

    # Define pixel
    out <- vector_build_grid(umap, N = 30, plot = plot)

    # Build network
    out <- vector_build_net(out, CUT = 1, plot = plot)

    # Calculate Quantile Polarization (QP) score
    out <- vector_build_value(out, pca, plot = plot)

    # Get pixel's QP score
    out <- vector_grid_value(out, plot = plot)

    # Find starting point
    out <- vector_auto_center(out, up_value = 0.9, plot = plot)

    # Infer vector
    out <- vector_draw_arrow(
      out,
      P = 0.9,
      plot = plot,
      colors = out$colors,
      plot.SUMMIT = plot
    )

    if (plot) {
      graphics::par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
      graphics::mtext("UMAP_1", side = 1)
      graphics::mtext("UMAP_2", side = 2)
    }

    return(out)
  }
)

#' @param filter Logical value, default set to \code{FALSE},
#' will filter cells though pseudotime.
#'
#' @method infer_vector default
#'
#' @rdname infer_vector
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' vector_result <- infer_vector(
#'   example_matrix,
#'   plot = FALSE
#' )
#' }
setMethod(
  f = "infer_vector",
  signature = signature(object = "matrix"),
  definition = function(
      object,
      filter = FALSE,
      plot = TRUE,
      ...) {
    dimension <- get_embedding(object)
    result <- infer_vector(dimension)
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
        pseudotime_vector = pseudotime
      )
    )

    return(object)
  }
)

#' @param reduction The reduction method.
#' @param dims The reduction to use.
#' 
#' @method infer_vector Seurat
#'
#' @rdname infer_vector
setMethod(
  f = "infer_vector",
  signature = signature(object = "Seurat"),
  definition = function(
      object,
      filter = FALSE,
      plot = TRUE,
      reduction = "umap",
      dims = 2,
      ...) {
    dimension <- get_embedding(
      object,
      reduction = reduction,
      dims = dims
    )
    result <- infer_vector(dimension)
    selected <- which(result$P.PS != "NA")
    pseudotime <- result$P.PS[selected]
    pseudotime <- normalization(pseudotime)
    object <- object[, selected]
    object@meta.data$pseudotime_vector <- pseudotime

    return(object)
  }
)
