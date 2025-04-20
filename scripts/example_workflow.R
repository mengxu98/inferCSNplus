

object <- initiate_object(example_matrix)
object <- inferCSN(object)

plot_gene_rank(object)






perturb_tfs <- c(
  "g1" = 0
)

object <- predictPerturbation(
  object = object,
  perturb_tfs = perturb_tfs,
  use_weight = TRUE,
  n_iter = 3,
  scale_factor = 0.1,
  verbose = TRUE
)

object <- embedPerturbation(
  object = object,
  reduction_method = "umap",
  dims = 2,
  scale = TRUE,
  cores = 1,
  seed = 42
)

p1 <- plotPerturbation(
  object = object,
  label_points = FALSE,
  point_size = 0.5,
  alpha = 0.5,
  show_trajectory = FALSE
)
p1

# 交替展示方式（默认）
p2 <- plotPerturbationTrajectory(
  object = object,
  genes = names(perturb_tfs),
  show_heatmap = TRUE,
  n_top_genes = 10,
  heatmap_style = "combined"
)
p2$heatmap
p2$trajectory

# 左右分开展示方式
p2_sep <- plotPerturbationTrajectory(
  object = object,
  genes = names(perturb_tfs),
  show_heatmap = TRUE,
  n_top_genes = 10,
  heatmap_style = "separate"
)
p2_sep$heatmap
p2_sep$trajectory

# 保存热图
pdf("heatmap_combined.pdf", width = 10, height = 12)
print(p2$heatmap)
dev.off()

pdf("heatmap_separate.pdf", width = 15, height = 12)  # 加宽以适应左右布局
print(p2_sep$heatmap)
dev.off()

pert_results <- list(
  original = object@perturbation@original,
  perturbed = object@perturbation@perturbed,
  delta = object@perturbation@perturbed - object@perturbation@original,
  trajectory = object@perturbation@params$reduction$trajectory_vectors
)

ggsave("perturbation_plot.pdf", p1, width = 8, height = 6)
ggsave("trajectory_plot.pdf", p2$trajectory, width = 10, height = 6)
if (!is.null(p2$heatmap)) {
  pdf("expression_heatmap.pdf", width = 8, height = 10)
  print(p2$heatmap)
  dev.off()
}
