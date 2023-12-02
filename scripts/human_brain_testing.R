library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
options(stringsAsFactors = FALSE)
options(timeout = 1000)

data_input_path <- "DIRECTNET/human_Brain_atac-seq/"
data_output_path <- "DIRECTNET/human_Brain_atac-seq_output/"

data <- readRDS(
  file = paste0(
    data_input_path,
    "GSE147672_200324_Brain_scATAC_SummarizedExperiment_forGEO"
  )
)

counts <- SummarizedExperiment::assay(data, "counts")
# Filter Chrs to ensure that the remaining Chrs expression level is greater than 1000
counts <- counts[which(rowSums(counts) > 1000), ]
meta <- as.data.frame(SummarizedExperiment::colData(data))

# import low dimension data and cluster information
umap_info <- readxl::read_excel(paste0(data_input_path,
                                       "41588_2020_721_MOESM5_ESM.xlsx"),
                                skip = 21)

cell_coord <- cbind(umap_info$UMAP_dim1, umap_info$UMAP_dim2)
rownames(cell_coord) <- paste0("scATAC_CTRL_",
                               umap_info$Region, "_",
                               umap_info$Donor_ID, "#",
                               umap_info$`10x_SingleCell_Barcode`, "-1")

colnames(cell_coord) <- c("UMAP_1", "UMAP_2")
index <- match(colnames(counts), rownames(cell_coord))
identical(colnames(counts), rownames(cell_coord)[index])
cell_coord <- cell_coord[index, ]
identical(colnames(counts), rownames(cell_coord))
group0 <- umap_info$Cluster[index]

# based on original nature genetics paper: https://www.nature.com/articles/s41588-020-00721-x
group <- c(
  "Striatal inhibitory 1",
  "Oligodendrocytes 2",
  "Striatal inhibitory 2",
  "Oligodendrocytes 1",
  "unclassified",
  "Oligodendrocytes 4",
  "Isocortical exitatory",
  "OPCs 2",
  "unclassified",
  "Oligodendrocytes 5",
  "Microglia", "unclassified",
  "Oligodendrocytes 3",
  "OPCs 1", "unclassified",
  "Hippocampal exitatory 2",
  "Isocortical astrocytes",
  "Striatal astrocytes",
  "Hippocampal exitatory 1",
  "Isocortical astrocytes",
  "unclassified", "Isocortical inhibitory",
  "unclassified", "Nigral OPCs"
)
group_new <- plyr::mapvalues(group0, from = unique(group0), to = group)

index <- which(group_new != "unclassified")
counts <- counts[, index]
meta <- meta[index, ]
meta$celltype <- group_new[index]

genome_info <- read.table(file = paste0(data_input_path,
                                        "hg38.promoter.regions.txt"))
names(genome_info) <- c("Chrom", "Starts", "Ends", "genes")
genes <- lapply(genome_info$genes, function(x) strsplit(x, "[|]")[[1]][1])
genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
genes <- unlist(genes)
genome_info$genes <- genes
unik <- !duplicated(genes) # filter out different transcript
genome_info <- genome_info[unik, ]

markers <- read.table(file = paste0(data_input_path,
                                    "Brain_markers.txt"),
                      sep = "\t")

markers <- data.frame(gene = markers$gene, group = markers$cluster)
clusters_unique <- unique(markers$group)
marker_list <- list()
for (i in seq_along(clusters_unique)) {
  marker1 <- markers[markers$group == clusters_unique[i], ]
  marker_list[[i]] <- as.character(marker1$gene)
}

# New ###############################
markers$annotation <- plyr::mapvalues(
  markers$group,
  from = unique(markers$group),
  to = group[which(group != "unclassified")]
)
# End ###############################

markers_all <- unique(unlist(marker_list))
focus_markers <- markers_all

# Part I: Set up a Seurat object
combined <- CreateSeuratObject(counts = counts,
                               meta.data = meta,
                               assay = "ATAC")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 50)
combined <- RunSVD(combined)
DepthCor(combined)

combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)
Idents(combined) <- combined$celltype
combined@reductions$umap@cell.embeddings <- cell_coord

# save(combined, file = paste0(data_output_path, "combined.rds"))
# load(paste0(data_output_path, "combined.rds"))

# Part II: Infer CREs using DIRECT-NET
combined <- Run_DIRECT_NET(
  combined,
  peakcalling = FALSE,
  k_neigh = 50,
  atacbinary = TRUE,
  max_overlap = 0.5,
  size_factor_normalize = TRUE,
  genome.info = genome_info,
  focus_markers = "FEV"
)

direct.net_result <- Misc(combined, slot = "direct.net")
direct.net_result <- as.data.frame(do.call(cbind, direct.net_result))

# We have run DIRECT-NET on all markers, load the results for downstream network analysis
# save(direct.net_result, file = paste0(data_output_path, "Brain_direct.net.RData"))
# load(paste0(data_output_path, "Brain_direct.net.RData"))

# check the function type name
direct.net_result$function_type <- gsub(
  "HF", "HC", direct.net_result$function_type
)
direct.net_result$function_type <- gsub(
  "Rest", "MC", direct.net_result$function_type
)
direct.net_result$function_type <- gsub(
  "LF", "LC", direct.net_result$function_type
)

# Part III: Visualize the links
temp <- tempfile()
download.file(
  "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz",
  temp
)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
Plot_connections(direct.net_result, gene_anno, marker = "FEV", cutoff = 0.2)

# Part IV: Construct regulatory networks
# identify differential accessible peaks, here we focus on Striatal inhibitory 1 and Striatal inhibitory 2
### identify differential accessible peaks (DA)
DefaultAssay(combined) <- "ATAC"

# focused_markers <- markers[which(markers$group %in% c("Striatal inhibitory 1", "Striatal inhibitory 2")), , drop = FALSE]
focused_markers <- markers[which(markers$annotation %in% c("Striatal inhibitory 1", "Striatal inhibitory 2")), , drop = FALSE]
# groups <- unique(focused_markers$group)
groups <- unique(focused_markers$annotation)
da_peaks_list <- list()
for (i in seq_along(groups)) {
  print(i)
  da_peaks <- FindMarkers(
    object = combined,
    min.pct = 0.2,
    ident.1 = groups[i],
    group.by = "celltype",
    test.use = "LR",
    only.pos = TRUE
  )
  da_peaks_list[[i]] <- da_peaks
}

# Detect CRE-TF connections
# CRE-gene connections
CREs_Gene <- generate_CRE_Gene_links(
  direct.net_result,
  markers = focused_markers
)
# Find focused CREs which is overlapped with DA
Focused_CREs <- generate_CRE(
  L_G_record = CREs_Gene$distal,
  P_L_G_record = CREs_Gene$promoter,
  da_peaks_list
)

# detect TFs for distal CREs
library(BSgenome.Hsapiens.UCSC.hg38)
L_TF_record <- generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$distal,
  species = "Homo sapiens",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  markers = focused_markers
)
# detect TFs for Promoters
P_L_TF_record <- generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$promoter,
  species = "Homo sapiens",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  markers = focused_markers
)

# Output TF-gene connections
network_links <- generate_links_for_Cytoscape(
  L_G_record = CREs_Gene$distal,
  L_TF_record,
  P_L_G_record = CREs_Gene$promoter,
  P_L_TF_record,
  groups
)

# Ouput Node attribute
Node_attribute <- generate_node_for_Cytoscape(
  network_links,
  markers = focused_markers
)
