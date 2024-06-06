library(SingleCellExperiment)
library(slingshot)

#-----------------------------------step1
scobj <- seurat_object_male

scale.data <- scobj@assays$RNA$scale.data
scale.gene <- rownames(scale.data)

counts <- scobj@assays$RNA$counts
counts <- counts[scale.gene, ]

sim <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts)
)

#----------------------------------step2
umap <- scobj@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP-1", "UMAP-2")
SingleCellExperiment::reducedDims(sim) <- S4Vectors::SimpleList(UMAP = umap)
plot(umap, pch = 16, asp = 1)

meta <- scobj@meta.data
colData(sim)$orig.ident <- meta$orig.ident
colData(sim)$celltype <- meta$CellType

#----------------------------------step3
sim <- slingshot::slingshot(
  sim,
  clusterLabels = "celltype",
  reducedDim = "UMAP"
)
colnames(colData(sim))

#-------------------------------step4
colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks = 100)]
plotcol[is.na(plotcol)] <- "lightgrey"
plotcol
sim@colData$plotcol <- plotcol
plot(reducedDims(sim)$UMAP, col = plotcol, pch = 16, asp = 1)
lines(
  slingshot::SlingshotDataSet(sim),
  lwd = 2,
  col = RColorBrewer::brewer.pal(9, "Set1")
)
legend(
  "right",
  legend = paste0("lineage", 1:6),
  col = unique(RColorBrewer::brewer.pal(6, "Set1")),
  pch = 16
)

#-----------------------------------step5 tradeSeq
# Fit negative binomial model
counts <- sim@assays@data$counts
crv <- slingshot::SlingshotDataSet(sim)
pseudotime <- slingshot::slingPseudotime(crv, na = FALSE)
cell_weights <- slingshot::slingCurveWeights(crv)

set.seed(1)
icMat <- tradeSeq::evaluateK(
  counts = counts,
  sds = crv,
  k = 3:10, # no more than 12
  nGenes = 500,
  parallel = TRUE,
  verbose = TRUE
)

# fit negative binomial GAM
# 2k cells ~13 min
counts <- counts[intersect(TFs, rownames(counts)), ]

sce <- tradeSeq::fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  nknots = 14,
  parallel = TRUE,
  verbose = TRUE
)

 # Exploring the correlation between gene expression and pseudo timing
asso_res <- tradeSeq::associationTest(sce)
head(asso_res)
# Searching for genes with the highest correlation with start and end points
start_res <- tradeSeq::startVsEndTest(sce)
head(start_res)
# Order and select the genes with the strongest correlation (waldStat) and visualize them
o_start <- order(start_res$waldStat, decreasing = TRUE)
sig_genes <- names(sce)[o_start[1]]
sig_genes <- "NR2E3"
tradeSeq::plotSmoothers(sce, counts, gene = sig_genes)

endRes <- tradeSeq::diffEndTest(sce,pairwise=T)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
sigGene <- "NR2E3"
tradeSeq::plotSmoothers(sce, counts, sigGene)
tradeSeq::plotGeneCount(crv, counts, gene = sigGene)

#Lineage1
tmp1 <- rownames(start_res [start_res$logFClineage1 < 0.01,])
tmp2 <- rownames(endRes[endRes$logFC1_2<0.01,])
intersect(tmp1, tmp2)

coldata <- data.frame(
  celltype = sim@colData$celltype,
  plotcol = sim@colData$plotcol
)
rownames(coldata) <- colnames(sim)

filter_coldata <- coldata[colnames(sce), ]

filter_coldata$Pseudotime <- sce$crv$pseudotime.Lineage1

top_genes <- names(sce)[o_start]
matrix_top_genes <- sce@assays@data$counts[top_genes, ]
matrix_top_genes <- log10(matrix_top_genes + 1)
hm_dyn(
  matrix_top_genes,
  dyn_tfs_male,
  topX = 50,
  # cRow = TRUE,
  # filename = "fig03_male.pdf",
  width = 4.5,
  limits = c(-10, 20),
  # toScale = T,
  height = 4.5
)


top6 <- names(sce)[o_start[1:14]]
top6_exp <- sce@assays@data$counts[top6, ]
top6_exp <- t(log2(top6_exp + 1))


plt_data <- cbind(filter_coldata, top6_exp)
colnames(plt_data)

library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
mycolors <- getPalette(length(unique(top6)))

plt_list <- list()
for (gene in top6) {
  print(gene)
  p <- ggscatter(
    data = plt_data,
    x = "Pseudotime",
    y = gene,
    color = "celltype",
    size = 0.6
  ) +
    geom_smooth(se = FALSE, color = "orange") +
    theme_bw() +
    scale_color_manual(values = mycolors) +
    theme(legend.position = "none")

  plt_list[[gene]] <- p
}

patchwork::wrap_plots(plt_list)
