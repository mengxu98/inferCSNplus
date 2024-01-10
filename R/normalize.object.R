#' @title normalize.seurat.object
#'
#' @param object Seurat object
#'
#' @return Normalized seurat object
#' @export
#'
normalize.seurat.object <- function(object) {
  if (!is(object, "Seurat")) {
    stop("Pleasure input an Seurat object")
  }

  if ("RNA" %in% names(object@assays)) {
    Seurat::DefaultAssay(object) <- "RNA"
    object <- Seurat::NormalizeData(
      object = object,
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )
    object <- Seurat::ScaleData(object = object)
    object <- Seurat::RunPCA(
      object,
      features = Seurat::VariableFeatures(object = object),
      pc.genes = object@var.genes,
      pcs.compute = 40,
      do.print = FALSE
    )
    object <- Seurat::RunUMAP(object, dims = 1:40)
    all_markers_list <- Seurat::FindAllMarkers(object)
    all_markers_list <- all_markers_list[all_markers_list$p_val_adj <= 0.05, ]
    all_markers_list <- all_markers_list[all_markers_list$avg_log2FC >= 1, ]
    Seurat::Misc(object, slot = "all_markers_list") <- all_markers_list
  }

  if ("ATAC" %in% names(object@assays)) {
    Seurat::DefaultAssay(object) <- "ATAC"
    object <- Signac::RunTFIDF(object)
    object <- Signac::FindTopFeatures(
      object,
      min.cutoff = "q0"
    )
    object <- Signac::RunSVD(object)
    object <- Seurat::RunUMAP(
      object,
      reduction = "lsi",
      dims = 2:50,
      reduction.name = "umap.atac",
      reduction.key = "atacUMAP_"
    )
  }

  if (all(c("RNA", "ATAC") %in% names(object@assays))) {
    object <- Seurat::FindMultiModalNeighbors(
      object,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:50, 2:50)
    )
    object <- Seurat::RunUMAP(
      object,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    )
  }
}
