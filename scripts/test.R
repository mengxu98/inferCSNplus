library(BSgenome.Hsapiens.UCSC.hg38)
library(inferCSN)
data(motifs)

seurat_object <- readRDS(
  "~/Study/research/_repositories/seurat_object_multiome_EN-neurons.RDS"
)
Seurat::Idents(seurat_object) <- "celltype"

seurat_object <- subset(
  seurat_object,
  idents = c("RG", "IPC", "EN", "EN-fetal-late", "EN-fetal-early")[1:2]
)


object <- initiate_object(
  seurat_object,
  filter_mode = "unfiltered", # "variable"
  filter_by = "celltype",
  peak_assay = "ATAC",
  celltype_by = "celltype"
)

object <- find_motifs(
  object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

object <- inferCSN(
  object,
  n_cores = 1
)

networks <- export_csn(
  object
)

rg_network <- export_csn(
  object,
  celltype = "RG",
  weight_cutoff = 0.1
)
head(rg_network)

DefaultNetwork(object)
GetNetwork(object)


tfs <- GetTFs(grn_summary)
length(tfs)

# Print inferred coefficients
coef(object)

# Find gene and regulatory modules
test_srt <- find_modules(object)

# Print modules
NetworkModules(object)
