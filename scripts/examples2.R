load("../stab.combined.RData")
meta_data <- h.combined@meta.data
matrix <- h.combined@assays$RNA@counts

sample_num <- 2000

meta_data_astro <- meta_data[meta_data$cluster == "Astro", ]
meta_data_astro_male <- meta_data_astro[meta_data_astro$Sex == "Male", ]
meta_data_astro_female <- meta_data_astro[meta_data_astro$Sex == "Female", ]
matrix_male <- matrix[, rownames(meta_data_astro_male)][, sample(10000, 2000)]
matrix_female <- matrix[, rownames(meta_data_astro_female)][, sample(2800, 2000)]

seurat_astro_male <- Seurat::CreateSeuratObject(
  counts = matrix_male,
  assay = "RNA",
  meta.data = meta_data_astro_male
)
seurat_astro_male <- FindVariableFeatures(seurat_astro_male)
seurat_astro_male <- NormalizeData(seurat_astro_male)
seurat_astro_male <- ScaleData(seurat_astro_male)
seurat_astro_male <- RunPCA(seurat_astro_male)
seurat_astro_male <- FindNeighbors(seurat_astro_male, dims = 1:20)
seurat_astro_male <- FindClusters(seurat_astro_male, resolution = 1)

seurat_astro_male <- RunUMAP(seurat_astro_male, dims = 1:20, return.model = T)
DimPlot(seurat_astro_male, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(seurat_astro_male, label = T, label.size = 4, reduction = "umap", group.by = "Sample")
DimPlot(seurat_astro_male, label = T, label.size = 3, reduction = "umap", group.by = "cluster") + NoLegend()
DimPlot(seurat_astro_male, label = T, label.size = 3, reduction = "umap", group.by = "cluster1") + NoLegend()
DimPlot(seurat_astro_male, label = T, label.size = 4, reduction = "umap", group.by = "Period")

seurat_astro_male <- Seurat::SCTransform(seurat_astro_male)
DefaultAssay(seurat_astro_male) <- "SCT"


# seurat_astro_male <- FindVariableFeatures(seurat_astro_male)
# seurat_astro_male <- NormalizeData(seurat_astro_male)
# seurat_astro_male <- ScaleData(seurat_astro_male)
# seurat_astro_male <- RunPCA(seurat_astro_male)
seurat_astro_male <- harmony::RunHarmony(seurat_astro_male, "Sample")
seurat_astro_male <- RunUMAP(seurat_astro_male, dims = 1:20, reduction = "harmony")
seurat_astro_male <- FindNeighbors(seurat_astro_male, dims = 1:20, reduction = "harmony")
seurat_astro_male <- FindClusters(seurat_astro_male, resolution = 1)
DimPlot(seurat_astro_male, label = T, label.size = 4, reduction = "umap", group.by = "Sample")
DimPlot(seurat_astro_male, label = T, label.size = 3, reduction = "umap", group.by = "cluster") + NoLegend()
DimPlot(seurat_astro_male, label = T, label.size = 3, reduction = "umap", group.by = "cluster1") + NoLegend()
DimPlot(seurat_astro_male, label = T, label.size = 4, reduction = "umap", group.by = "Period")



seurat_astro_male_list <- SplitObject(seurat_astro_male, split.by = "Sample")

seurat_astro_male_list <- lapply(
  X = seurat_astro_male_list, FUN = function(x) {
    if (ncol(x) < 100) {
      return()
    }
    x <- NormalizeData(x)
    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 2000
    )
  }
)
seurat_astro_male_list <- seurat_astro_male_list[!map_lgl(seurat_astro_male_list, is.null)]

features <- SelectIntegrationFeatures(object.list = seurat_astro_male_list)
seurat_astro_male_anchors <- FindIntegrationAnchors(object.list = seurat_astro_male_list, anchor.features = features)
seurat_astro_male_combined <- IntegrateData(anchorset = seurat_astro_male_anchors)
DefaultAssay(seurat_astro_male_combined) <- "integrated"
seurat_astro_male_combined <- ScaleData(seurat_astro_male_combined, verbose = FALSE)
seurat_astro_male_combined <- RunPCA(seurat_astro_male_combined, npcs = 30, verbose = FALSE)
seurat_astro_male_combined <- RunUMAP(seurat_astro_male_combined, reduction = "pca", dims = 1:30)
seurat_astro_male_combined <- FindNeighbors(seurat_astro_male_combined, reduction = "pca", dims = 1:30)
seurat_astro_male_combined <- FindClusters(seurat_astro_male_combined)

DimPlot(seurat_astro_male_combined, label = T, label.size = 4, reduction = "umap", group.by = "Sample")
DimPlot(seurat_astro_male_combined, label = T, label.size = 3, reduction = "umap", group.by = "cluster") + NoLegend()
DimPlot(seurat_astro_male_combined, label = T, label.size = 3, reduction = "umap", group.by = "cluster1") + NoLegend()
DimPlot(seurat_astro_male_combined, label = T, label.size = 4, reduction = "umap", group.by = "Period")


seurat_astro_male1 <- inferVECTOR(
  seurat_astro_male_combined
)


seurat_astro_female <- Seurat::CreateSeuratObject(
  counts = matrix_female,
  assay = "RNA",
  meta.data = meta_data_astro_female
)

seurat_astro_female <- FindVariableFeatures(seurat_astro_female)
seurat_astro_female <- NormalizeData(seurat_astro_female)
seurat_astro_female <- ScaleData(seurat_astro_female)
seurat_astro_female <- RunPCA(seurat_astro_female)
seurat_astro_female <- FindNeighbors(seurat_astro_female, dims = 1:50)
seurat_astro_female <- FindClusters(seurat_astro_female, resolution = 1)

seurat_astro_female <- RunUMAP(seurat_astro_female, dims = 1:50, return.model = T)
DimPlot(seurat_astro_female, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(seurat_astro_female, label = T, label.size = 3, reduction = "umap", group.by = "cluster") + NoLegend()
DimPlot(seurat_astro_female, label = T, label.size = 3, reduction = "umap", group.by = "cluster1") + NoLegend()
DimPlot(seurat_astro_female, label = T, label.size = 4, reduction = "umap", group.by = "Period")
DimPlot(seurat_astro_female, label = T, label.size = 4, reduction = "umap", group.by = "Sample")


seurat_astro_female <- Seurat::SCTransform(seurat_astro_female)
DefaultAssay(seurat_astro_female) <- "SCT"

seurat_astro_female <- harmony::RunHarmony(seurat_astro_female, "Sample")
DimPlot(seurat_astro_female, label = T, label.size = 4, reduction = "harmony", group.by = "Sample")

vector_result_female1 <- inferVECTOR(seurat_astro_female)
