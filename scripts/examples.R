object <- readRDS("../allen_m1c_2019_ssv4.rds")
tfs_list <- read.table("../regulators.txt")[, 1]

DefaultAssay(object) <- "RNA"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object)
object <- FindNeighbors(object, dims = 1:10)
object <- FindClusters(object, resolution = 0.5)

object <- RunUMAP(object, dims = 1:10, return.model = T)
DimPlot(object, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(object, label = T, label.size = 3, reduction = "umap", group.by = "subclass") + NoLegend()

object_male <- subset(object, donor_sex == "M")
object_male <- FindVariableFeatures(object_male)
object_male <- ScaleData(object_male)
object_male <- RunPCA(object_male)
object_male <- FindNeighbors(object_male, dims = 1:10)
object_male <- FindClusters(object_male, resolution = 0.5)

object_male <- RunUMAP(object_male, dims = 1:10, return.model = T)
DimPlot(object_male, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(object_male, label = T, label.size = 3, reduction = "umap", group.by = "subclass") + NoLegend()

object_male_new <- inferVECTOR(object_male)
DimPlot(object_male_new, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(object_male_new, label = T, label.size = 3, reduction = "umap", group.by = "subclass") + NoLegend()
DimPlot(object_male_new, label = T, label.size = 3, reduction = "umap", group.by = "region")
table(object_male_new$subclass)
table(object_male_new$region)
Seurat::Idents(object_male_new) <- "subclass"

k_neigh <- 50
atacbinary <- TRUE
max_overlap <- 0.8
reduction_name <- NULL
size_factor_normalize <- FALSE
genome_info <- NULL
verbose <- TRUE

agg_data <- aggregating.data(
  object_male_new,
  k_neigh = k_neigh,
  atacbinary = atacbinary,
  max_overlap = max_overlap,
  reduction_name = NULL,
  size_factor_normalize = size_factor_normalize,
  verbose = verbose
)
Seurat::Misc(object_male_new, slot = "aggregated_data") <- agg_data

object_male_new <- inferCSN(
  object_male_new,
  regulators = tfs_list,
  verbose = TRUE,
  cores = 1,
  aggregate = FALSE
)

object_female <- subset(object, donor_sex == "F")
object_female <- FindVariableFeatures(object_female)
object_female <- ScaleData(object_female)
object_female <- RunPCA(object_female)
object_female <- FindNeighbors(object_female, dims = 1:10)
object_female <- FindClusters(object_female, resolution = 0.5)

object <- RunUMAP(object_female, dims = 1:10, return.model = T)
DimPlot(object_female, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(object_female, label = T, label.size = 3, reduction = "umap", group.by = "subclass") + NoLegend()

object_female_new <- inferVECTOR(object_female)
DimPlot(object_female_new, label = T, label.size = 3, reduction = "umap") + NoLegend()
DimPlot(object_female_new, label = T, label.size = 3, reduction = "umap", group.by = "subclass") + NoLegend()
table(object_female_new$subclass)

object_female_new <- inferCSN(
  object_female_new,
  regulators = tfs_list,
  verbose = TRUE,
  cores = 1,
  aggregate = FALSE
)
