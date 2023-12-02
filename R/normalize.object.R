normalize.object <- function(object) {
  if (class(object) != "Seurat") {
    stop("Pleasure input an Seurat object")
  }
  if (all(c("RNA", "ATAC") %in% names(object@assays))) {
    # RNA analysis
    DefaultAssay(combined) <- "RNA"
    combined <- NormalizeData(object = combined,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

    #####
    # Detection of variable genes across the single cells
    combined <- FindVariableFeatures(object = combined)
    features.rna <- markers$gene
    combined@assays$RNA@var.features <- features.rna #Seurat 4 error
    # combined@assays[["RNA"]]@meta.data[["var.features.rank"]] <- features.rna # Seurat 5
    combined <- Seurat::FindMarkers(object = combined)
    #####

    ##### New
    all.genes <- rownames(combined)
    combined <- ScaleData(combined, features = all.genes)
    combined <- RunPCA(combined, features = VariableFeatures(object = combined))
    #####
    combined <- ScaleData(object = combined)
    combined <- RunPCA(combined, pc.genes = combined@var.genes, pcs.compute = 40, do.print = FALSE)
    combined <- RunUMAP(combined, dims = 1:40)

    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(combined) <- "ATAC"
    combined <- Signac::RunTFIDF(combined)
    #> Performing TF-IDF normalization
    combined <- Signac::FindTopFeatures(combined, min.cutoff = 'q0')
    combined <- Signac::RunSVD(combined)
    combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    #
    combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


    Idents(combined) <- combined$celltype

  }

}
