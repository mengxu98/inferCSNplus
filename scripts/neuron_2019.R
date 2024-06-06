# https://doi.org/10.1016/j.neuron.2019.06.011
devtools::document()

source("scripts/functions.R")
library(Seurat)

load("~/HAR/brain_data/neuron-sc_dev_cortex_geschwind/raw_counts_mat.rdata")
meta_data <- read.csv(
  "~/HAR/brain_data/neuron-sc_dev_cortex_geschwind/cell_metadata.csv",
  row.names = 1
)

seurat_object <- Seurat::CreateSeuratObject(
  raw_counts_mat,
  assay = "RNA",
  meta.data = meta_data
)

save(seurat_object, file = "../seurat_object_neuron.RData")

# new mission
load(file = "../seurat_object_neuron.RData")

seurat_object <- seurat_object[, !is.na(seurat_object$Cluster)]

seurat_object <- seurat_object[, which(seurat_object$Subcluster != c("ExN_5"))]
seurat_object <- seurat_object[, which(seurat_object$Subcluster != c("ExM_7"))]
seurat_object <- seurat_object[, which(seurat_object$Subcluster != c("ExN_6"))]
seurat_object <- seurat_object[, which(seurat_object$Subcluster != c("ExM-U_5"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("OPC"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("Mic"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("End"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("InCGE"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("InMGE"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("Per"))]
seurat_object <- seurat_object[, which(seurat_object$Cluster != c("PgG2M"))]
seurat_object$Cluster[which(seurat_object$Cluster == "PgG2M")] <- "PG"
seurat_object$Cluster[which(seurat_object$Cluster == "PgS")] <- "PG"
seurat_object$Cluster[which(seurat_object$Cluster == "vRG")] <- "RG"
seurat_object$Cluster[which(seurat_object$Cluster == "oRG")] <- "RG"

# test <- select_Seurat(seurat_object, type = "unselect")
# seurat_object <- seurat_object[, test]
seurat_object <- seurat_function(
  seurat_object,
  group_by = "Cluster",
  dims_num = 30,
  resolution = 1
)

Seurat::DimPlot(
  seurat_object,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "Cluster"
)

Seurat::DimPlot(
  seurat_object,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "Subcluster"
)

seurat_object <- get.pseudotime(
  seurat_object,
  cluster_by = "Cluster",
  start_cluster = "RG"
)
seurat_object <- inferVECTOR(seurat_object)

p1 <- Seurat::DimPlot(
  seurat_object,
  cols = color_list,
  label = TRUE,
  label.size = 3,
  reduction = "umap",
  group.by = "Cluster"
) + Seurat::NoLegend() +
  Seurat::FeaturePlot(seurat_object, "pseudotime", cols = c("white", "#3366ff")) +
  Seurat::FeaturePlot(seurat_object, "Lineage1", cols = c("white", "#cc0033"))
p1
ggplot2::ggsave(filename = "../neuron_2019/fig01.pdf", p1, height = 3, width = 10)

Lineage <- "Lineage1"
meta_data <- seurat_object@meta.data
meta_data <- meta_data[which(!is.na(meta_data[[Lineage]])), ]
count <- seurat_object@assays$RNA$data[, rownames(meta_data)]
# count <- seurat_object@assays$RNA$scale.data[, rownames(meta_data)]

seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[Lineage]])]

seurat_object <- Seurat::FindVariableFeatures(seurat_object, nfeatures = 3000)
genes <- Seurat::VariableFeatures(seurat_object)

cor(
  meta_data$Lineage1,
  meta_data$pseudotime
)

# https://doi.org/10.1016/j.neuron.2019.01.027
# the introduction of this paper maybe could explane this trend
density_points(
  meta_data = meta_data,
  group_column = "Cluster",
  pseudotime_column = Lineage,
  min_cells = 50,
  plot = TRUE
)
ggplot2::ggsave(filename = "../neuron_2019/fig02.pdf", height = 4, width = 6)

dynamic_object_list <- dynamic.genes_new(
  object = count,
  meta_data = meta_data,
  group_column = "Cluster",
  pseudotime_column = Lineage,
  min_cells = 50,
  p_value = 0.01,
  cores = 10
)

epoch_assignments_list <- lapply(
  dynamic_object_list, function(x) {
    assign_network(
      matrix = count,
      dynamic_object = x
    )
  }
)

TFs <- read.table("scripts/regulators.txt")[, 1]
TFs_bd <- read.table("../BrainSpan/TFs_brain-development.txt", header = TRUE, sep = "\t")
TFs_bd <- TFs_bd$gene.symbol
TFs <- intersect(TFs, TFs_bd)

dynamic_genes <- dynamic.genes(
  object = count,
  pseudotime = meta_data$Lineage1,
  p_value = 0.01,
  cores = 10
)

coding_genes <- read.table("scripts/CCDS-genes.txt", header = TRUE)[, 1]
dynamic_genes <- intersect(dynamic_genes, coding_genes)
dynamic_genes2 <- intersect(dynamic_genes, genes)


save(list=ls(), file = "../neuron_2019/current_environment_data.Rdata")
load(file = "../neuron_2019/current_environment_data.Rdata")

grn_list <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    dyn_genes <- unique(x$mean_expression$gene)
    cells <- y[[2]]$cells$cells
    inferCSN(
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dyn_genes,
      cores = 6,
      verbose = TRUE
    )
  }
)
grn_list_0.5k_0.5r <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    dyn_genes <- unique(x$mean_expression$gene)
    cells <- y[[2]]$cells$cells
    inferCSN(
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dyn_genes,
      k_folds = 0.5,
      r_threshold = 0.5,
      cores = 6,
      verbose = TRUE
    )
  }
)

grn_list1 <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    cells <- y[[2]]$cells$cells
    inferCSN(
      # t(as.matrix(count[genes, cells])),
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dynamic_genes,
      cores = 6,
      verbose = TRUE
    )
  }
)

grn_list1_0.5k_0.5r <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    cells <- y[[2]]$cells$cells
    inferCSN(
      # t(as.matrix(count[genes, cells])),
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dynamic_genes,
      k_folds = 0.5,
      r_threshold = 0.5,
      cores = 6,
      verbose = TRUE
    )
  }
)

grn_list2 <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    cells <- y[[2]]$cells$cells
    inferCSN(
      # t(as.matrix(count[genes, cells])),
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dynamic_genes2,
      cores = 6,
      verbose = TRUE
    )
  }
)

grn_list2_0.5k_0.5r <- purrr::map2(
  epoch_assignments_list,
  dynamic_object_list,
  function(x, y) {
    cells <- y[[2]]$cells$cells
    inferCSN(
      # t(as.matrix(count[genes, cells])),
      t(as.matrix(count[, cells])),
      regulators = TFs,
      targets = dynamic_genes2,
      k_folds = 0.5,
      r_threshold = 0.5,
      cores = 6,
      verbose = TRUE
    )
  }
)
r2()

save(list=ls(), file = "../neuron_2019/current_environment_data.Rdata")
load(file = "../neuron_2019/current_environment_data.Rdata")

grn_list_top <- grn_list2_0.5k_0.5r
grn_list_top <- purrr::map(
  grn_list1_0.5k_0.5r,
  .f = function(x) {
    nsd <- 1
    thr <- mean(abs(x$weight)) + (nsd * sd(abs(x$weight)))
    message(thr)
    x[abs(x$weight) > thr, ]
  }
)
colnames_order <- gtools::mixedsort(names(grn_list_top))
lda_data <- construct_lda_data(grn_list_top, celltypes_order = colnames_order)
lda_result <- lda_analysis(
  lda_data, k = 10, path = "../neuron_2019/lda/", ntop = 30 ,
  binary = TRUE,
  binary_threshold = 0
) # , binary = TRUE

colnames_order_all <- c(
  "network_1..1_network_1..2",
  "network_1..2_network_2..2",
  "network_2..2_network_2..3",
  "network_2..3_network_3..3",
  "network_3..3_network_3..4",
  "network_3..4_network_4..4",
  "network_4..4_network_4..5",
  "network_4..5_network_5..5"
)
lda_result_sel <- lda_result[, colnames_order_all]
lda_result_list <- purrr::map(
  seq_len(ncol(lda_result_sel)),
  .f = function(x) {
    lda_result_sel[, x]
  }
)
tfs_intersect <- purrr::reduce(lda_result_list, intersect)

networks_order <- c(
  "network_1..1",
  "network_1..2",
  "network_2..2",
  "network_2..3",
  "network_3..3",
  "network_3..4",
  "network_4..4",
  "network_4..5",
  "network_5..5"
)


grn_list_top_tfs <- purrr::map_dfr(
  networks_order,
  .f = function(x) {
    network <- grn_list_top[x][[1]]
    network <- network_format(
      network,
      regulators = tfs_intersect,
      abs_weight = FALSE
    )
    network$celltype <- x
    as.data.frame(network)
  }
)

plot_dynamic_networks(
  grn_list_top_tfs,
  celltypes_order = networks_order,
  figure_save = T,
  figure_name = "../neuron_2019/network.pdf",
  width = 12,
  height = 3.5,
  nrow = 1,
  ntop = 50,
  # layout = "kamadakawai",
  theme_type = "theme_facet"
)


grn_list_top_tfs_all <- purrr::map_dfr(
  colnames_order,
  .f = function(x) {
    network <- grn_list_top[x][[1]]
    network <- network_format(
      network,
      regulators = tfs_intersect,
      abs_weight = FALSE
    )
    network$celltype <- x
    as.data.frame(network)
  }
)
plot_dynamic_networks(
  grn_list_top_tfs_all,
  celltypes_order = colnames_order,
  figure_save = T,
  figure_name = "../neuron_2019/network.mp4",
  width = 12,
  height = 3.5,
  nrow = 1,
  ntop = 50,
  # layout = "kamadakawai",
  theme_type = "theme_facet",
  plot_type = "animate"
)


network_data <- purrr::map_dfr(
  colnames_order[1:6],
  .f = function(x) {
    network <- grn_list2[x][[1]][1:500, ]
    network$celltype <- x
    as.data.frame(network)
  }
)

plot_dynamic_networks(
  network_data,
  celltypes_order = colnames_order[1:6],
  figure_save = T,
  width = 12,
  height = 10
)


# save(grn_list, file = "../grn_list.Rdata")
# load("../grn_list.Rdata")

transition_networks_list <- extract_genes_transition(
  network_list = grn_list,
  top_regulators_num = 1,
  top_targets_num = 10000,
  common_regulators = TRUE,
  common_targets = TRUE
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0(
        "../enrich_results_KEGG3/1-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"
      ),
      figure_save = TRUE,
      plot = TRUE,
      enrich_method = "KEGG",
      method = "enrichR"
    )
  }
}


for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results_CellType2/1-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
      figure_save = TRUE,
      plot = TRUE,
      enrich_method = "CellType",
      method = "enrichR"
    )
  }
}

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0(
        "../enrich_results_GO2/1-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"
      ),
      figure_save = TRUE,
      plot = TRUE
    )
  }
}














transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 1,
  top_targets_num = 50
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/1-regulator_50-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
      plot = TRUE,
      figure_save = TRUE
    )
  }
}


transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 1000,
  top_targets_num = 500000
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/all-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
      plot = TRUE,
      figure_save = TRUE
    )
  }
}


transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 100,
  top_targets_num = 50
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/all-regulator_50-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
      plot = TRUE,
      figure_save = TRUE
    )
  }
}



transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 100,
  top_targets_num = 50,
  common_regulators = TRUE,
  common_targets = FALSE
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/all-regulator_50-targets_common_regulators&not_common_targets/trans_net_", i, "_", namess, ".png"),
      plot = TRUE,
      figure_save = TRUE
    )
  }
}



transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 100,
  top_targets_num = 50,
  common_regulators = FALSE,
  common_targets = FALSE
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/all-regulator_50-targets_not_common_regulators&not_common_targets/trans_net_", i, "_", namess, ".png"),
      plot = TRUE,
      figure_save = TRUE
    )
  }
}



Idents(seurat_object) <- seurat_object$CellType
pbmc.markers.all <- Seurat::FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

pbmc.markers.all_sub <- dplyr::filter(pbmc.markers.all, gene == genes_list)
selectTFsNums <- length(rankgene)
selectTFsNums <- 13

pbmc.markers.all_sub <- pbmc.markers.all[1:selectTFsNums, ]
pbmc.markers.all_sub$gene <- rankgene[1:selectTFsNums]


pbmc.markers <- pbmc.markers.all_sub %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = selectTFsNums, wt = avg_log2FC)

# check
head(pbmc.markers)


# write.csv(pbmc.markers, "pbmc.markers.csv")


# prepare data from seurat object
st.data <- prepareDataFromscRNA(
  object = pbmc,
  diffData = pbmc.markers,
  showAverage = TRUE
)





# Plot a heatmap of the dynamic TFs
dgenes <- purrr::map_dfr(
  epoch_assignments_list, function(x) {
    x$mean_expression
  }
)
dgenes <- unique(dgenes$gene)

tfstoplot <- intersect(dgenes, TFs)
dyn_tfs <- list(
  meta_data = meta_data,
  genes = tfstoplot
)

meta_data$cells <- rownames(meta_data)
# First, smooth expression for a cleaner plot
count_smooth <- expression_ksmooth(
  count,
  meta_data = meta_data,
  pseudotime_column = "Lineage1",
  bandwith = 0.03
)

matrix_hm <- seurat_object@assays$RNA$counts
matrix_hm <- pre_pseudotime_matrix(matrix_hm)

hm_dyn(
  # count,
  pt.matrix,
  # count_smooth,
  # matrix_hm,
  dyn_tfs,
  topX = 50,
  # filename = "fig03.png",
  width = 4.5,
  limits = c(-10, 20),
  toScale = T,
  height = 4.5
)

# Plot a heatmap of all dynamic TFs and target genes
# dyngenes <- xdyn
# dyngenes$genes <- dyngenes$genes[names(dynTFs$genes) %in% dgenes]
# hm_dyn(count_smooth, dyngenes, topX=100)

grn_list_1000 <- lapply(
  grn_list, function(x) {
    x[1:50000, ]
  }
)

tfs_list <- purrr::map(
  grn_list_1000, function(x) {
    x[, 1]
  }
)


tfs_list <- purrr::map(
  grn_list, function(x) {
    unique(x[, 1])
  }
)
tfs_intersect <- purrr::reduce(
  tfs_list,
  .f = intersect
)
tfs_intersect <- TFs

plot_dynamic_network(
  grn_list,
  tfs_intersect,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig04.png", height = 10, width = 10)


plot_dynamic_network(
  list(mesoderm_network = grn_data),
  TFs,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig05.png", height = 10, width = 10)

plot_detail_network(
  grn_list,
  regulators = tfs_intersect,
  top_edges = 100,
  only_TFs = TRUE,
  network_order = c("network_1..1", "network_1..2", "network_2..2"),
  communities = NULL,
  compute_betweenness = TRUE
)
# ggplot2::ggsave(filename = "fig06.png", height = 10, width = 10)

gene_rank <- compute_pagerank(grn_list)

plot_top_features(
  grn_list,
  gene_ranks = gene_rank,
  regulators = tfs_intersect,
  targets = tfs_intersect,
  # top_regulators_num = 10,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_top_features(
  grn_list,
  gene_ranks = gene_rank,
  # regulators = tfs_intersect,
  # targets = tfs_intersect[1:100],
  top_regulators_num = 20,
  top_targets_num = 10,
  method = "page_rank",
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_top_features(
  grn_list,
  regulators = NULL,
  targets = tfs_intersect[1:10],
  regulators_num = 1,
  targets_num = NULL,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_targets_with_top_regulators_detail(
  grn_list,
  epochs_list = epoch_assignments_list,
  targets = tfs_intersect[1:10],
  regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_targets_and_regulators(
  grn_list,
  epochs = epoch_assignments_list,
  targets = tfs_intersect[1:10],
  regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_targets_with_top_regulators(
  grn_list,
  gene_ranks = gene_rank,
  targets = tfs_intersect[1:10],
  top_regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)


cell_metadata <- seurat_object@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(seurat_object)
)
rownames(gene_annotation) <- rownames(seurat_object)
matrix <- Seurat::GetAssayData(seurat_object, layer = "count", assay = "RNA")
cds <- new_cell_data_set(matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

save(cds, file = "../cds_neuron.RData")

# new mission
load(file = "../cds_neuron.RData")

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

plot_cells(cds, color_cells_by = "Cluster", group_label_size = 4)
## Step 5: Learn a graph
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "Cluster", group_label_size = 4)

## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5
)



plot_cells(cds)
