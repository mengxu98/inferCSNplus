#' @title annotation.celltype
#'
#' @param object seurat object
#' @param method
#'
#' @return seurat object
#' @export
#'
celltype.annotation <- function(
    # https://github.com/Teichlab/celltypist
    # celltypist$models$download_models(force_update = T) # First run needs to download the trained model data.
    # source_python('celltypist_model_download.py')
    # rownames(object@assays$RNA@counts)
    object,
    method = "celltypist",
    reference_path = NULL) {
  message(paste0("[", Sys.time(), "] -----: Run '", method, "' for cell type annotation."))
  if (method == "celltypist") {
    pandas <- reticulate::import("pandas")
    numpy <- reticulate::import("numpy")
    scanpy <- reticulate::import("scanpy")
    celltypist <- reticulate::import("celltypist")

    # matrix <- Matrix::t(Matrix::as.matrix(object[["RNA"]]@counts))
    matrix <- Matrix::t(Matrix::as.matrix(object@assays$RNA$counts))
    meta_data <- object@meta.data

    adata <- scanpy$AnnData(
      X = numpy$array(matrix),
      obs = pandas$DataFrame(meta_data),
      var = pandas$DataFrame(data.frame(
        gene = rownames(object@assays$RNA$counts),
        row.names = rownames(object@assays$RNA$counts)
      ))
    )

    # adata <- scanpy$AnnData(
    #   X = numpy$array(matrix),
    #   obs = pandas$DataFrame(meta_data),
    #   var = pandas$DataFrame(
    #     data.frame(
    #       gene = colnames(matrix),
    #       row.names = colnames(matrix)
    #     )
    #   )
    # )

    model <- celltypist$models$Model$load(model = reference_path) # Developing_Human_Brain.pkl
    model$cell_types
    scanpy$pp$normalize_total(adata, target_sum = 1e4)
    scanpy$pp$log1p(adata)
    predictions <- celltypist$annotate(adata, model = reference_path, majority_voting = FALSE)
    object <- AddMetaData(object, predictions$predicted_labels)
    object$celltype <- object$predicted_labels
    return(object)
  } else if (method == "SingleR") {
    hpca.se <- SingleR::HumanPrimaryCellAtlasData() # Bug
    # load(annotation.celltype)
    matrix <- object@assays$RNA@counts
    singler_obj <- SingleR::SingleR(test = matrix, ref = hpca.se, labels = hpca.se$label.main)
    object$celltype <- singler_obj$labels
    return(object)
  } else if (!(method %in% c("celltypist", "SingleR"))) {
    message("[", Sys.time(), "] -----: Please choose method in: 'celltypist' or 'SingleR!'")
  }
}
