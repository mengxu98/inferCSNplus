devtools::document()

source("scripts/functions.R")
library(Seurat)

data_path <- "~/HAR/brain_data/GSE162170/"
data_path <- "../../"

RNA_SCE <- readRDS(paste0(data_path, "Multiome_RNA_SCE.RDS"))
ATAC_SCE <- readRDS(paste0(data_path, "Multiome_ATAC_SCE.RDS"))

# ATAC_GeneActivity <- readRDS(paste0(data_path, "Multiome_ATAC_GeneActivity.RDS"))

genes_data <- as.data.frame(RNA_SCE@rowRanges@elementMetadata@listData)
genes_data_protein_coding <- dplyr::filter(genes_data, gene_type == "protein_coding")

rna_count <- RNA_SCE@assays@data$counts
rna_count <- rna_count[genes_data_protein_coding$gene_id_short, ]
rownames(rna_count) <- genes_data_protein_coding$gene_name
rna_count <- rna_count[!duplicated(rownames(rna_count)), ]

meta_data <- read.table(
  paste0(data_path, "multiome_cell_metadata.txt"),
  sep = "\t",
  header = TRUE
)
cluster_data <- read.table(
  paste0(data_path, "multiome_cluster_names..txt"),
  sep = "\t",
  header = TRUE
)
meta_data_rna <- as.data.frame(RNA_SCE@colData)
meta_data_atac <- as.data.frame(ATAC_SCE@colData)

cluster_data_rna <- dplyr::filter(cluster_data, Assay == "Multiome RNA")
cluster_data_atac <- dplyr::filter(cluster_data, Assay == "Multiome ATAC")

colnames(cluster_data_rna) <- c("Assay", "seurat_clusters", "Cluster.Name")

meta_data_new <- purrr::map2_dfr(
  cluster_data_rna$seurat_clusters,
  cluster_data_rna$Cluster.Name,
  .f = function(x, y) {
    b <- subset(meta_data_rna, seurat_clusters == x)
    b$celltype <- y
    b
  }
)

seurat_object <- Seurat::CreateSeuratObject(
  rna_count,
  assay = "RNA",
  meta.data = meta_data_new
)

atac_matrix <- ATAC_SCE@assays@data$counts
ranges <- ATAC_SCE@rowRanges

annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"
genome.name <- "hg38"
GenomeInfoDb::genome(annotations) <- genome.name

chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_matrix,
  ranges = ranges,
  sep = c(":", "-"),
  genome = genome.name,
  min.cells = 0,
  annotation = annotations
)

# combine
seurat_object[["ATAC"]] <- chrom_assay
save(seurat_object, file = "../seurat_object.RData")

# new mission
load(file = "../seurat_object.RData")

seurat_object <- subset(
  seurat_object,
  DF_classification == "Singlet"
)

seurat_object <- seurat_object[, which(seurat_object$celltype != c("SP"))]
seurat_object <- seurat_object[, which(seurat_object$celltype != c("IN3"))]
seurat_object <- seurat_object[, which(seurat_object$celltype != c("Cyc. Prog."))]
seurat_object <- seurat_object[, which(seurat_object$celltype != c("mGPC/OPC"))]

seurat_object <- seurat_function(
  seurat_object,
  group_by = "celltype",
  dims_num = 35,
  resolution = 0.5
)

seurat_object <- inferVECTOR(seurat_object)
seurat_object <- get.pseudotime(
  seurat_object,
  cluster_by = "celltype",
  start_cluster = "RG"
)
seurat_object <- seurat_object[, !is.na(seurat_object$Lineage1)]

DimPlot(
  seurat_object,
  # cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "celltype"
) + NoLegend()

cor(
  seurat_object$Lineage1,
  seurat_object$pseudotime
)

# Select variable features
seurat_object <- FindVariableFeatures(
  seurat_object,
  assay = "RNA"
)

genes <- VariableFeatures(seurat_object, assay = "RNA")

genes <- dynamic.genes(
  seurat_object,
  pseudotime_column = "Lineage1",
  cluster_column = "celltype",
  cores = 6
)

# running
# BiocManager::install("TFBSTools")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("motifmatchr")
library(BSgenome.Hsapiens.UCSC.hg38)

# Get motif data
data(motifs)

# Initiate GRN object and select candidate regions
seurat_object <- initiate_object(
  seurat_object,
  peak_assay = "ATAC"
)

# Scan candidate regions for TF binding motifs
seurat_object <- find_motifs(
  seurat_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer gene regulatory network
seurat_object <- inferCSN(
  seurat_object,
  cores = 12,
  genes = genes[1:100],
  method = "srm",
  verbose = 2,
  cores = 1
)

save(seurat_object, file = "../seurat_object_grn.RData")
load(file = "../seurat_object_grn.RData")

seurat_object_list <- initiate_object_list.CSNObject(
  seurat_object,
  cluster_column = "celltype",
  peak_assay = "ATAC",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
seurat_object_list <- inferCSN(seurat_object_list, verbose = TRUE)

# Print inferred coefficients
net <- coef(seurat_object)
net

# Find gene and regulatory modules
seurat_object <- find_modules(
  seurat_object,
  p_thresh = 0.05,
  rsq_thresh = 0.3,
  nvar_thresh = 10,
  min_genes_per_module = 5
)

plot_gof(seurat_object, point_size = 2)
plot_module_metrics(seurat_object)
seurat_object <- get_network_graph(
  seurat_object
)
plot_network_graph(seurat_object)
plot_network_graph(seurat_object, layout = "fr")


genes_num <- 100
method <- "srm"
bench::mark(
  inferCSN2(
    seurat_object,
    genes = genes[1:genes_num],
    method = method
  ),
  inferCSN2(
    seurat_object,
    genes = genes[1:genes_num],
    method = method,
    parallel = FALSE,
    verbose = FALSE
  ),
  inferCSN2(
    seurat_object,
    genes = genes[1:genes_num],
    method = method,
    parallel = TRUE
  ),
  inferCSN(
    seurat_object,
    cores = 1,
    genes = genes[1:genes_num],
    method = method,
    verbose = 2
  ),
  inferCSN(
    seurat_object,
    cores = 6,
    genes = genes[1:genes_num],
    method = method,
    verbose = 2
  ),
  inferCSN(
    seurat_object,
    cores = 12,
    genes = genes[1:genes_num],
    method = method,
    verbose = 2
  ),
  memory = FALSE,
  check = FALSE
)
