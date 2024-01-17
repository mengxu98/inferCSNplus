library(Seurat)
# DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
seurat_object <- CreateSeuratObject(counts = DATA, project = "seurat_object3k", min.cells = 0, min.features = 0)
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object),npcs = 150)
seurat_object <- RunUMAP(seurat_object, dims = 1:50)


DimPlot(seurat_object, reduction = "umap")
saveRDS(seurat_object, file='seurat_object.RDS')

seurat_object_sub <- seurat_object
seurat_object_sub <- subset(seurat_object, subclass == "Astro")
DimPlot(seurat_object_sub, reduction = "umap")

VEC = seurat_object_sub@reductions$umap@cell.embeddings
rownames(VEC) = colnames(seurat_object_sub)
PCA = seurat_object_sub@reductions$pca@cell.embeddings
rownames(PCA) = colnames(seurat_object_sub)

# Remove quantile-based colinearity among PCs (new feature in VECTOR 0.0.3):
PCA=vector.rankPCA(PCA)

# Define pixel
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Quantile Polarization (QP) score
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's QP score
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find starting point
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)

# OUT$P.PS : Peseudotime Score (PS) of each cell

