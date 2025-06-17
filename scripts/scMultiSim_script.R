library(scMultiSim)

data(GRN_params_100)
grn <- GRN_params_100
grn$regulated.gene <- paste0("g_", grn$regulated.gene)
grn$regulator.gene <- paste0("g_", grn$regulator.gene)

sim_results <- sim_true_counts(
  list(
    rand.seed = 0,
    GRN = grn,
    num.cells = 500,
    num.cifs = 50,
    cif.sigma = 0.5,
    tree = Phyla1(),
    diff.cif.fraction = 0.8,
    do.velocity = TRUE,
    dynamic.GRN = list(
      cell.per.step = 3,
      num.changing.edges = 2,
      weight.mean = 0,
      weight.sd = 4
    )
  )
)
names(sim_results)

grn_gd <- grn
colnames(grn_gd) <- c("target", "regulator", "weight")
grn_gd <- grn_gd[, c("target", "regulator", "weight")]
grn_gd <- network_format(grn_gd, abs_weight = F)

counts <- sim_results$counts
counts_t <- t(counts)
network <- inferCSN(
  counts_t,
  penalty = "L0L2",
  cores = 4,
  r_squared_threshold = 0
)
calculate_metrics(
  network,
  grn_gd[, 1:2],
  return_plot = TRUE
)

cell_meta <- results$cell_meta

network_genie3 <- GENIE3::GENIE3(
  counts,
  verbose = TRUE,
  nCores = 4
)
network_genie3 <- GENIE3::getLinkList(
  network_genie3
)
colnames(network_genie3) <- c("target", "regulator", "weight")
network_genie3 <- network_format(network_genie3, abs_weight = F)
calculate_metrics(
  network_genie3,
  grn_gd[, 1:2],
  return_plot = TRUE
)

plot_tsne(log2(results$counts + 1),
  results$cell_meta$pop,
  legend = "pop", plot.name = "True RNA Counts Tsne"
)

plot_tsne(log2(results$atacseq_data + 1),
  results$cell_meta$pop,
  legend = "pop", plot.name = "True ATAC-seq Tsne"
)

plot_rna_velocity(results, arrow.length = 2)

plot_gene_module_cor_heatmap(results)

add_expr_noise(
  results,
  # options go here
  alpha_mean = 1e4
)
names(results)
results$cell_time

counts_obs <- results$counts_obs
counts_obs1 <- results$counts

plot_grn(grn)
