source("scripts/functions.R")
library(Seurat)

# load("data/hsTFs.rda")
# TFs <- intersect(rownames(count), hsTFs)
TFs <- read.table("scripts/regulators.txt")[, 1]

seurat_list <- list()
seurat_list[[1]] <- readRDS("../brain_data/GSE192774/GSE192773_Seurat_Microglia_RNA.RDS")
seurat_list[[2]] <- readRDS("../brain_data/GSE192774/GSE192773_Seurat_Excitatory_RNA.RDS")
seurat_list[[3]] <- readRDS("../brain_data/GSE192774/GSE192773_Seurat_OPCOli_RNA.RDS")
seurat_list[[4]] <- readRDS("../brain_data/GSE192774/GSE192773_Seurat_Astrocyte_RNA.RDS")
seurat_list[[5]] <- readRDS("../brain_data/GSE192774/GSE192773_Seurat_Inhibitory_RNA.RDS")

seurat_object <- merge(
  seurat_list[[1]],
  y = seurat_list[2:5],
  project = "brain"
)
seurat_object <- seurat_list[[2]]

########
seurat_object <- readRDS("/Users/mx/Downloads/GSE192773_Seurat_Excitatory_RNA.RDS")

seurat_object <- subset(seurat_object, Species == "human")

seurat_object$Sex <- "F"
seurat_object$Sex[seurat_object$Sample == "h2"] <- "M"
seurat_object$Sex[seurat_object$Sample == "h3"] <- "M"

seurat_object_male <- subset(seurat_object, Sex == "M")
seurat_object_male <- seurat_function(seurat_object_male)

seurat_object_female <- subset(seurat_object, Sex == "F")
seurat_object_female <- seurat_function(seurat_object_female)


p1 <- DimPlot(
  seurat_object_male,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap"
) + NoLegend()

p2 <- DimPlot(
  seurat_object_male,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "Sample"
) + NoLegend()

p3 <- DimPlot(
  seurat_object_male,
  cols = color_list[-c(1:4)],
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "CellType"
) + NoLegend()


p <- p1 + p2 + p3
p
ggplot2::ggsave(p, filename = "fig01_male.pdf", width = 12, height = 4)

# pdf("fig02_male.pdf", width = 9, height = 6)
# seurat_object_male <- inferVECTOR(object = seurat_object_male)
# dev.off()

seurat_object_male <- get.pseudotime(
  object = seurat_object_male,
  cluster_by = "CellType",
  assay = "integrated",
  slot = "data"
)

meta_data_male <- seurat_object_male@meta.data
meta_data_male <- meta_data_male[which(!is.na(meta_data_male$Lineage1)), ]
count_male <- seurat_object_male@assays$integrated@data[, rownames(meta_data_male)]

dynamic_object_list <- findDynGenes_new(
  count_male,
  meta_data,
  group_column = "CellType",
  pseudotime_column = "Lineage1"
)

epoch_assignments_male_list <- lapply(
  dynamic_object_list, function(x){
    assign_epochs_new(
      matrix = count_male,
      dynamic_object = x
    )
  }
)

grn_list <- list()
for (i in 1:7) {
  x <- epoch_assignments_male_list[[i]]
  y <- dynamic_object_list[[i]]
  grn_list[[i]] <- inferCSN(
    t(count_male[, y[[2]]$cells$cells]),
    regulators = TFs,
    targets = x$mean_expression$gene,
    cores = 6,
    verbose = TRUE
  )
}

grn_list <- purrr::map2(
  epoch_assignments_male_list,dynamic_object_list, function(x, y) {
    inferCSN(
      t(count_male[, y[[2]]$cells$cells]),
      regulators = TFs,
      targets = x$mean_expression$gene,
      cores = 6,
      verbose = TRUE
    )
  }
)

grn_list <- lapply(
  epoch_assignments_male_list, function(x) {
    inferCSN(
      t(count_male),
      regulators = TFs,
      targets = x$mean_expression$gene,
      cores = 6,
      verbose = TRUE
    )
  }
)


# First, smooth expression for a cleaner plot
cells <- meta_data
count_smooth_male <- expression_ksmooth(count_male, cells, bandwith = 0.03)

# Plot a heatmap of the dynamic TFs
dgenes_male <- purrr::map_dfr(
  epoch_assignments_male_list, function(x) {
    x$mean_expression
  }
)
dgenes_male <- dgenes_male$gene

tfstoplot_male <- intersect(dgenes_male, TFs)
dyn_tfs_male <- xdyn_male
dyn_tfs_male$genes <- dyn_tfs_male$genes[
  names(dyn_tfs_male$genes) %in% tfstoplot_male
]

hm_dyn(
  count_smooth_male,
  dyn_tfs_male,
  topX = 50,
  # filename = "fig03_male.pdf",
  width = 4.5,
  height = 4.5
)

# Plot a heatmap of all dynamic TFs and target genes
# dyngenes <- xdyn
# dyngenes$genes <- dyngenes$genes[names(dynTFs$genes) %in% dgenes]
# hm_dyn(count_smooth, dyngenes, topX=100)

plot_dynamic_network_new(
  grn_list,
  TFs,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig04_male.pdf", height = 10, width = 10)


plot_dynamic_network(
  list(mesoderm_network = grn_data_male),
  TFs,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig05_male.pdf", height = 10, width = 10)

plot_dynnet_detail(
  dynamic_grn_male,
  tfs = TFs,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3"),
  communities = NULL,
  compute_betweenness = TRUE
)
ggplot2::ggsave(filename = "fig06_male.pdf", height = 10, width = 10)

plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  TFs,
  numTopTFs = 10,
  only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  tfstoplot_male,
  numTopTFs = 10, only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

# We can specify additional parameters including the number of top TFs and targets:
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  TFs,
  numTopTFs = 3,
  numTargets = 5,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)









dynamic_grn_male_list <- purrr::map2(
  grn_list, epoch_assignments_male_list, function(x, y) {
    epochGRN(
      x,
      y
    )
  }
)

dynamic_grn_male <- epochGRN(
  grn_data_male,
  epoch_assignments_male
)
dynamic_grn_male <- as.data.frame(dynamic_grn_male)
gene_rank_male <- compute_pagerank(
  dynamic_grn_male,
  weight_column = "weighted_score"
)


grn_list <- list()
for (cluster in unique(meta_data_male$CellType)) {
  # Find dynamically expressed genes
  xdyn_male <- findDynGenes(
    count_male,
    meta_data = meta_data_male,
    cluster_by = cluster,
    group_column = "CellType",
    pseudotime_column = "Lineage1"
  )

  dgenes_male <- na.omit(names(xdyn_male$genes)[xdyn_male$genes < 0.01])
  count_male_sub <- count_male[dgenes_male, xdyn_male$cells$cells]

  if (dim(count_male_sub)[1] > 2 && dim(count_male_sub)[1] > 100) {
    grn_list[[cluster]] <- inferCSN(
      t(count_male_sub),
      regulators = TFs,
      cores = 6,
      verbose = TRUE
    )
  } else {
    grn_list[[cluster]] <- NULL
  }
}

dynamic_object <- xdyn_male_all

xdyn_male_all <- findDynGenes(
  count_male,
  meta_data = meta_data_male,
  group_column = "CellType",
  pseudotime_column = "Lineage1"
)

xdyn_male_all <- define_epochs(
  xdyn_male_all,
  count_male,
  method = "pseudotime",
  num_epochs = 5
)

epoch_assignments_male <- assign_epochs(
  matrix = count_male,
  dynamic_object = xdyn_male_all,
  forceGenes = FALSE
)

dgenes_male_all <- na.omit(names(xdyn_male_all$genes)[xdyn_male_all$genes < 0.01])
# Reconstruct and perform crossweighting
grn_data_male <- reconstructGRN(
  count_male,
  TFs,
  dgenes_male_all,
  method = "pearson",
  zThresh = 3
)

# # so slow
# grn_data_male <- crossweight(
#   grn_data_male,
#   count_male,
#   xdyn_male,
#   filter_thresh = 0.5
# )

dynamic_grn_male <- epochGRN(
  grn_data_male,
  epoch_assignments_male
)
dynamic_grn_male <- as.data.frame(dynamic_grn_male)
gene_rank_male <- compute_pagerank(
  dynamic_grn_male,
  weight_column = "corr"
)

# First, smooth expression for a cleaner plot
ccells_male <- xdyn_male$cells
count_smooth_male <- expression_ksmooth(count_male, ccells_male, bandwith = 0.03)

# Plot a heatmap of the dynamic TFs
tfstoplot_male <- intersect(dgenes_male, TFs)
dyn_tfs_male <- xdyn_male
dyn_tfs_male$genes <- dyn_tfs_male$genes[
  names(dyn_tfs_male$genes) %in% tfstoplot_male
]

save(
  dynamic_grn_male,
  count_smooth_male,
  dyn_tfs_male, TFs,
  grn_data_male,
  gene_rank_male,
  tfstoplot_male,
  file = "plot_data_male.Rdata"
)

hm_dyn(
  count_smooth_male,
  dyn_tfs_male,
  topX = 50,
  filename = "fig03_male.pdf",
  width = 4.5,
  height = 4.5
)

# Plot a heatmap of all dynamic TFs and target genes
# dyngenes <- xdyn
# dyngenes$genes <- dyngenes$genes[names(dynTFs$genes) %in% dgenes]
# hm_dyn(count_smooth, dyngenes, topX=100)

plot_dynamic_network(
  dynamic_grn_male,
  TFs,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
)
ggplot2::ggsave(filename = "fig04_male.pdf", height = 10, width = 10)


plot_dynamic_network(
  list(mesoderm_network = grn_data_male),
  TFs,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig05_male.pdf", height = 10, width = 10)

plot_dynnet_detail(
  dynamic_grn_male,
  tfs = TFs,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3"),
  communities = NULL,
  compute_betweenness = TRUE
)
ggplot2::ggsave(filename = "fig06_male.pdf", height = 10, width = 10)

plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  TFs,
  numTopTFs = 10,
  only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  tfstoplot_male,
  numTopTFs = 10, only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

# We can specify additional parameters including the number of top TFs and targets:
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  TFs,
  numTopTFs = 3,
  numTargets = 5,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)






seurat_object_female <- seurat_function(seurat_object_female)

p1 <- DimPlot(
  seurat_object_female,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap"
) + NoLegend()

p2 <- DimPlot(
  seurat_object_female,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "Sample"
) + NoLegend()

p3 <- DimPlot(
  seurat_object_female,
  cols = color_list[-c(1:4)],
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "CellType"
) + NoLegend()

p <- p1 + p2 + p3
ggplot2::ggsave(p, filename = "fig01_female.pdf", width = 12, height = 4)

pdf("fig02_female.pdf", width = 9, height = 6)
seurat_object_female_vector <- inferVECTOR(object = seurat_object_female)
dev.off()

count_female <- as.matrix(seurat_object_female_vector@assays$RNA$data)
meta_data_female <- seurat_object_female_vector@meta.data

# Find dynamically expressed genes
xdyn_female <- findDynGenes(
  count_female,
  meta_data_female,
  group_column = "CellType",
  pseudotime_column = "pseudotime"
)

dgenes_female <- na.omit(names(xdyn_female$genes)[xdyn_female$genes < 0.01])

# dgenes <- dynamic.genes(object = Microglia_RNA)

# Reconstruct and perform crossweighting
grn_data_female <- reconstructGRN(
  count_female,
  TFs,
  dgenes_female,
  method = "pearson",
  zThresh = 3
)

# so slow
grn_data_female <- crossweight(
  grn_data_female,
  count_female,
  xdyn_female,
  filter_thresh = 0.5
)

xdyn_female <- define_epochs(
  xdyn_female,
  count_female,
  method = "pseudotime",
  num_epochs = 3
)
epoch_assignments_female <- assign_epochs(count_female, xdyn_female)

dynamic_grn_female <- epochGRN(grn_data_female, epoch_assignments_female)
dynamic_grn_female <- as.data.frame(dynamic_grn_female)
gene_rank_female <- compute_pagerank(
  dynamic_grn_female,
  weight_column = "weighted_score"
)

# First, smooth expression for a cleaner plot
ccells_female <- xdyn_female$cells
count_smooth_female <- expression_ksmooth(count_female, ccells_female, bandwith = 0.03)

# Plot a heatmap of the dynamic TFs
tfstoplot_female <- intersect(dgenes_female, TFs)
dyn_tfs_female <- xdyn_female
dyn_tfs_female$genes <- dyn_tfs_female$genes[
  names(dyn_tfs_female$genes) %in% tfstoplot_female
]

save(
  dynamic_grn_female,
  count_smooth_female,
  dyn_tfs_female, TFs,
  grn_data_female,
  gene_rank_female,
  tfstoplot_female,
  file = "plot_data_female.Rdata"
)

hm_dyn(
  count_smooth_female,
  dyn_tfs_female,
  topX = 50,
  filename = "fig03_female.pdf",
  width = 4.5,
  height = 4.5
)

genes <- names(sort(dyn_tfs_female$genes, decreasing = FALSE))[1:50]
plot_dynnet_detail(
  dynamic_grn_female,
  tfs = genes,
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3"),
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

genes <- names(sort(dyn_tfs_male$genes, decreasing = FALSE))[1:50]
plot_dynnet_detail(
  dynamic_grn_male,
  tfs = genes,
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3"),
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

plot_dynnet_detail(
  dynamic_grn_female,
  tfs = c(
    # "GLIS3",
    # "E2F3",
    # "ETS1"
    # "FOXO1",
    "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
  # order = c("epoch1..epoch2", "epoch2..epoch3")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)


plot_dynnet_detail(
  dynamic_grn_male,
  tfs = c(
    # "GLIS3",
    # "E2F3",
    # "ETS1"
    # "FOXO1",
    "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
  # order = c("epoch1..epoch2", "epoch2..epoch3")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

plot_dynnet_detail(
  dynamic_grn_female,
  tfs = c(
    "SOX5",
    "JUND",
    "TBL1X",
    "TCF4",
    "ENO1"
    # "NR2F1",
    # "DBP",
    # "RORA",
    # "NR1D1"
    # "FOXO1",
    # "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
  # order = c("epoch1..epoch2", "epoch2..epoch3")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

plot_dynnet_detail(
  dynamic_grn_male,
  tfs = c(
    "TCF4",
    "ZEB1",
    "BARX2",
    "PKNOX2",
    "SREBF2"
    # "FOXO1",
    # "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
  # order = c("epoch1..epoch2", "epoch2..epoch3")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)


plot_dynnet_detail(
  dynamic_grn_male,
  tfs = TFs[1:50],
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
  # order = c("epoch1..epoch2", "epoch2..epoch3")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

plot_top_regulators(
  dynamic_grn_female,
  gene_rank_female,
  c(
    # "GLIS3",
    # "E2F3",
    # "ETS1"
    "BMP7",
    "FOXO1",
    "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  numTopTFs = 2,
  only_TFs = FALSE,
  # order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  c(
    # "GLIS3",
    # "E2F3",
    # "ETS1"
    "FOXO1",
    "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  numTopTFs = 2,
  only_TFs = FALSE,
  # order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)
plot_top_regulators(
  dynamic_grn_male,
  gene_rank_male,
  c(
    "SOX5",
    "JUND",
    "TBL1X"
    # "FOXO1",
    # "FOXP2"
    # "SOX5",
    # "JUND",
    # "TBL1X"
  ),
  numTopTFs = 2,
  only_TFs = FALSE,
  # order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)
# Plot a heatmap of all dynamic TFs and target genes
# dyngenes <- xdyn
# dyngenes$genes <- dyngenes$genes[names(dynTFs$genes) %in% dgenes]
# hm_dyn(count_smooth, dyngenes, topX=100)

# plot_dynamic_network(
#   dynamic_grn_female,
#   TFs,
#   only_TFs = TRUE,
#   order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3")
# )
#
# plot_dynamic_network(
#   list(mesoderm_network = grn_data_female),
#   TFs,
#   only_TFs = TRUE
# )

plot_dynnet_detail(
  dynamic_grn_female,
  tfs = TFs,
  only_TFs = TRUE,
  communities = NULL,
  compute_betweenness = TRUE,
  # order = c("epoch1..epoch1", "epoch2..epoch2", "epoch3..epoch3"),
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

plot_top_regulators(
  dynamic_grn_female,
  gene_rank_female,
  TFs,
  numTopTFs = 10,
  only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)

plot_top_regulators(
  dynamic_grn_female,
  gene_rank_female,
  tfstoplot_female,
  numTopTFs = 10,
  only_TFs = FALSE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2", "epoch2..epoch3", "epoch3..epoch3")
)

# We can specify additional parameters including the number of top TFs and targets:
plot_top_regulators(
  dynamic_grn_female,
  gene_rank_female,
  TFs,
  numTopTFs = 3,
  numTargets = 5,
  only_TFs = TRUE,
  order = c("epoch1..epoch1", "epoch1..epoch2", "epoch2..epoch2")
)
