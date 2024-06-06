load("/Users/mx/Downloads/GitHub/MTG_object_sex.Rdata")
pal <- unique(MTG_object_female$cluster_color)[-1][-2][-2] # paletteer::paletteer_d("ggsci::nrc_npg")[1:16]

MTG_object_male <- inferVECTOR(MTG_object_male)

p1 <- DimPlot(
  MTG_object_male,
  cols= pal,
  label = T,
  label.size = 3,
  reduction = "umap",
  # group.by = "cluster_label"
) + NoLegend()

p2 <- DimPlot(
  MTG_object_male1,
  cols= pal,
  label = T,
  label.size = 3.5,
  reduction = "umap",
  group.by = "cluster_label"
)

p1 + p2


MTG_object_female <- inferVECTOR(MTG_object_female)

p3 <- DimPlot(
  MTG_object_female,
  cols= pal,
  label = T,
  label.size = 3,
  reduction = "umap",
  # group.by = "cluster_label"
) + NoLegend()

p4 <- DimPlot(
  MTG_object_female1,
  cols= pal,
  label = T,
  label.size = 3.5,
  reduction = "umap",
  group.by = "cluster_label"
)

p3 + p4

seurat_object_filter <- function(
    object,
    group,
    num = 50) {
  object_list <- Seurat::SplitObject(object, split.by = "cluster")
  for (i in 1:length(object_list)) {
    if (ncol(object_list[[i]]) < num) {
      object_list[[i]] <- NA
    }
  }
  object_list <- object_list[!map_lgl(object_list, is.na)]
  object <- merge(
    x = object_list[[1]],
    y = object_list[2:length(object_list)]
  )

  return(object)
}

tf_list <- read.table('../regulators.txt')[, 1]

Seurat::Idents(MTG_object_female) <- "cluster_label"
MTG_object_female$cluster <- Seurat::Idents(MTG_object_female)
# MTG_object_female1 <- seurat_object_filter(MTG_object_female1, "cluster")
MTG_object_female <- inferCSN(MTG_object_female, regulators = tf_list, verbose = TRUE)

object <- subset(MTG_object_female, cluster == "Exc L3-4 RORB CARM1P1")
object <- inferCSN(
  object,
  regulators = tf_list,
  aggregate = FALSE,
  verbose = TRUE
)
