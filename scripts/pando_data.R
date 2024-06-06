
RNA_ATAC <- readRDS("~/HAR/brain_data/Pando-data/seurat_objects/RNA_ATAC_metacells_srt.rds")

DimPlot(
  RNA_ATAC,
  # cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "pseudotime_ranks"
) + NoLegend()
