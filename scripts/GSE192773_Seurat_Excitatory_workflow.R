source("scripts/functions.R")
library(Seurat)

TFs <- read.table("scripts/regulators.txt")[, 1]

seurat_object <- readRDS(
  "../brain_data/GSE192774/GSE192773_Seurat_Excitatory_RNA.RDS"
)
seurat_object <- subset(seurat_object, Species == "human")

seurat_object$Sex <- "F"
seurat_object$Sex[seurat_object$Sample == "h2"] <- "M"
seurat_object$Sex[seurat_object$Sample == "h3"] <- "M"

seurat_object_male <- subset(
  seurat_object,
  Sex == "M"
)
seurat_object_male <- seurat_function(
  seurat_object_male,
  re_create_object = TRUE,
  group_by = "Sample",
  integrate = TRUE
)

seurat_object_female <- subset(
  seurat_object,
  Sex == "F"
)
seurat_object_female <- seurat_function(
  seurat_object_female,
  re_create_object = TRUE,
  group_by = "Sample",
  integrate = TRUE
)


seurat_object_male <- get.pseudotime(
  object = seurat_object_male,
  cluster_by = "CellType",
  assay = "RNA",
  slot = "data"
)

meta_data_male <- seurat_object_male@meta.data
meta_data_male <- meta_data_male[which(!is.na(meta_data_male$Lineage1)), ]
count_male <- seurat_object_male@assays$RNA$data[, rownames(meta_data_male)]

density_points(
  meta_data = meta_data_male,
  group_column = "CellType",
  pseudotime_column = "Lineage1",
  min_cells = 500,
  plot = TRUE
)

dynamic_object_list <- dynamic.genes_new(
  object = count_male,
  meta_data = meta_data_male,
  group_column = "CellType",
  pseudotime_column = "Lineage1",
  min_cells = 500,
  p_value = 0.01
)

epoch_assignments_male_list <- lapply(
  dynamic_object_list, function(x) {
    assign_network(
      matrix = count_male,
      dynamic_object = x
    )
  }
)

save(
  count_male,
  seurat_object_male,
  seurat_object_female,
  meta_data_male,
  meta_data,
  epoch_assignments_male_list,
  dynamic_object_list,
  TFs,
  file = "grn_data.Rdata"
)
load("grn_data.Rdata")

grn_list <- purrr::map2(
  epoch_assignments_male_list, dynamic_object_list, function(x, y) {
    genes <- unique(x$mean_expression$gene)
    cells <- y[[2]]$cells$cells
    inferCSN(
      t(count_male[genes, cells]),
      regulators = TFs,
      targets = genes,
      cores = 10,
      verbose = TRUE
    )
  }
)

save(grn_list, file = "grn_list.Rdata")
load("grn_list.Rdata")

transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 10,
  top_targets_num = 100,
  common_regulators = TRUE,
  common_targets = FALSE
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets[1:50],
      # filename = paste0("../enrich_results/1-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
      # figure_save = TRUE,
      plot = TRUE,
      enrich_method = "KEGG",
      method = "enrichR"
    )
  }
}


transition_networks_list <- extract_genes_transition(
  grn_list,
  top_regulators_num = 1,
  top_targets_num = 10000
)

for (i in 1:length(transition_networks_list)) {
  networks <- transition_networks_list[[i]]
  networks_names <- names(networks)
  for (j in 1:length(networks_names)) {
    namess <- networks_names[j]
    enrich_functions(
      networks[[namess]]$targets,
      filename = paste0("../enrich_results/1-regulator_all-targets_common_regulators&common_targets/trans_net_", i, "_", namess, ".png"),
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



Idents(seurat_object_male) <- seurat_object_male$CellType
pbmc.markers.all <- Seurat::FindAllMarkers(
  seurat_object_male,
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
dgenes_male <- purrr::map_dfr(
  epoch_assignments_male_list, function(x) {
    x$mean_expression
  }
)
dgenes_male <- unique(dgenes_male$gene)

tfstoplot_male <- intersect(dgenes_male, TFs)
dyn_tfs_male <- list(
  meta_data = meta_data,
  genes = tfstoplot_male
)

meta_data_male$cells <- rownames(meta_data_male)
# First, smooth expression for a cleaner plot
count_smooth_male <- expression_ksmooth(
  count_male,
  meta_data = meta_data_male,
  pseudotime_column = "Lineage1",
  bandwith = 0.03
)

matrix_hm <- seurat_object_male@assays$RNA$counts
matrix_hm <- pre_pseudotime_matrix(matrix_hm)

hm_dyn(
  # count_male,
  pt.matrix,
  # count_smooth_male,
  # matrix_hm,
  dyn_tfs_male,
  topX = 50,
  # filename = "fig03_male.png",
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
ggplot2::ggsave(filename = "fig04_male.png", height = 10, width = 10)


plot_dynamic_network(
  list(mesoderm_network = grn_data_male),
  TFs,
  only_TFs = TRUE
)
ggplot2::ggsave(filename = "fig05_male.png", height = 10, width = 10)

plot_detail_network(
  grn_list,
  regulators = tfs_intersect,
  top_edges = 100,
  only_TFs = TRUE,
  network_order = c("network_1..1", "network_1..2", "network_2..2"),
  communities = NULL,
  compute_betweenness = TRUE
)
# ggplot2::ggsave(filename = "fig06_male.png", height = 10, width = 10)

gene_rank_male <- compute_pagerank(grn_list)

plot_top_features(
  grn_list,
  gene_ranks = gene_rank_male,
  regulators = tfs_intersect,
  targets = tfs_intersect,
  # top_regulators_num = 10,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_top_features(
  grn_list,
  gene_ranks = gene_rank_male,
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
  epochs_list = epoch_assignments_male_list,
  targets = tfs_intersect[1:10],
  regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_targets_and_regulators(
  grn_list,
  epochs = epoch_assignments_male_list,
  targets = tfs_intersect[1:10],
  regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)

plot_targets_with_top_regulators(
  grn_list,
  gene_ranks = gene_rank_male,
  targets = tfs_intersect[1:10],
  top_regulators_num = 1,
  network_order = c("network_1..1", "network_1..2", "network_2..2")
)
