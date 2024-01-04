#' @title normalize.seurat.object
#'
#' @param object Seurat object
#'
#' @return Normalized seurat object
#' @export
#'
normalize.seurat.object <- function(object) {
  if (class(object) != "Seurat") {
    stop("Pleasure input an Seurat object")
  }

  if ("ATAC" %in% names(object@assays)) {
    DefaultAssay(object) <- "ATAC"
    object <- Signac::RunTFIDF(object)
    object <- Signac::FindTopFeatures(
      object,
      min.cutoff = "q0"
    )
    object <- Signac::RunSVD(object)
    object <- RunUMAP(
      object,
      reduction = 'lsi',
      dims = 2:50,
      reduction.name = "umap.atac",
      reduction.key = "atacUMAP_"
    )
  }

  if ("RNA" %in% names(object@assays)) {
    DefaultAssay(object) <- "RNA"
    object <- NormalizeData(
      object = object,
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )
    object <- ScaleData(object = object)
    object <- RunPCA(
      object,
      features = VariableFeatures(object = object),
      pc.genes = object@var.genes,
      pcs.compute = 40,
      do.print = FALSE
    )
    object <- RunUMAP(object, dims = 1:40)
    all_markers_list <- Seurat::FindAllMarkers(object)
    all_markers_list <- all_markers_list[all_markers_list$p_val_adj <= 0.05, ]
    all_markers_list <- all_markers_list[all_markers_list$avg_log2FC >= 1, ]
    Seurat::Misc(object, slot = "all_markers_list") <- all_markers_list
  }

  if (all(c("RNA", "ATAC") %in% names(object@assays))) {
    object <- FindMultiModalNeighbors(
      object,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:50, 2:50))
    object <- RunUMAP(
      object,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    )
  }
}
