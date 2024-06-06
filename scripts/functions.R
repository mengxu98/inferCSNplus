color_list <- c(
  "#3366cc", "#66ff66", "#3399ff", "#ffffcc", "#66ffff", "#66ccff",
  "#66cc66", "#ffcc66", "#66cccc", "#ffcccc", "#003366", "#ffcc00",
  "#ff9900", "#669966", "#ff9966", "#6699cc", "#ff99cc", "#6699ff",
  "#ff6600", "#666666", "#ff6666", "#6666cc", "#ff66cc", "#6666ff",
  "#ff3300", "#663366", "#ff3366", "#6633cc", "#ff33cc", "#6633ff",
  "#ff0000", "#660066", "#ff0066", "#6600cc", "#ff00cc", "#6600ff",
  "#000099", "#33ff66", "#ccff66", "#66ff99", "#ccffcc", "#33ffff",
  "#ffccff", "#33cc66", "#cccc66", "#66cc99", "#cccccc", "#33ccff",
  "#ff99ff", "#339966", "#cc9966", "#669999", "#cc99cc", "#66ffcc",
  "#ff66ff", "#336666", "#cc6666", "#666699", "#cc66cc", "#3366ff",
  "#ff33ff", "#333366", "#cc3366", "#663399", "#cc33cc", "#3333ff",
  "#ff00ff", "#330066", "#cc0066", "#660099", "#cc00cc", "#3300ff",
  "#99ff00", "#00ff66", "#99ff66", "#00ffcc", "#99ffcc", "#00ffff",
  "#99cc00", "#00cc66", "#99cc66", "#00cccc", "#99cccc", "#00ccff",
  "#999900", "#009966", "#999966", "#0099cc", "#9999cc", "#0099ff",
  "#996600", "#006666", "#996666", "#0066cc", "#9966cc", "#0066ff",
  "#993300", "#ffff00", "#993366", "#0033cc", "#9933cc", "#0033ff",
  "#990000", "#000066", "#990066", "#0000cc", "#9900cc", "#0000ff",
  "#66ff00", "#ffff99", "#33ffcc", "#ffff33", "#ccff00", "#66ff33",
  "#66cc00", "#ffcc99", "#33cccc", "#ffcc33", "#cccc00", "#66cc33",
  "#669900", "#ff9999", "#3399cc", "#ff9933", "#cc9900", "#669933",
  "#666600", "#ff6699", "#ffff66", "#ff6633", "#cc6600", "#666633",
  "#663300", "#ff3399", "#3333cc", "#ff3333", "#cc3300", "#663333",
  "#660000", "#ff0099", "#3300cc", "#ff0033", "#cc0000", "#660033",
  "#33ff00", "#ccff99", "#33ff99", "#ccff33", "#ccffff", "#33ff33",
  "#33cc00", "#cccc99", "#33cc99", "#cccc33", "#ccccff", "#33cc33",
  "#339900", "#cc9999", "#339999", "#cc9933", "#cc99ff", "#339933",
  "#336600", "#cc6699", "#336699", "#cc6633", "#cc66ff", "#336633",
  "#333300", "#cc3399", "#333399", "#cc3333", "#cc33ff", "#333333",
  "#330000", "#cc0099", "#330099", "#cc0033", "#cc00ff", "#330033",
  "#00ff00", "#99ff99", "#00ff99", "#99ff33", "#99ffff", "#00ff33",
  "#00cc00", "#99cc99", "#00cc99", "#99cc33", "#99ccff", "#00cc33",
  "#009900", "#999999", "#009999", "#999933", "#9999ff", "#009933",
  "#006600", "#996699", "#006699", "#996633", "#9966ff", "#006633",
  "#003300", "#993399", "#003399", "#993333", "#9933ff", "#003333",
  "#000000", "#990099", "#ffffff", "#990033", "#9900ff", "#000033"
)

seurat_function <- function(
    object,
    re_create_object = FALSE,
    dims_num = 50,
    integrate = FALSE,
    reduction = "umap",
    resolution = 1,
    group_by = NULL,
    plot = TRUE) {
  dims_num <- c(1:dims_num)
  if (re_create_object) {
    object <- Seurat::CreateSeuratObject(
      counts = object@assays$RNA$counts,
      assay = "RNA",
      meta.data = object@meta.data
    )
  }

  Seurat::DefaultAssay(object) <- "RNA"
  object <- Seurat::FindVariableFeatures(object)
  object <- Seurat::NormalizeData(object)
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object)
  object <- Seurat::FindNeighbors(object, dims = dims_num)
  object <- Seurat::FindClusters(object, resolution = resolution)
  object <- Seurat::RunUMAP(object, dims = dims_num)

  if (plot) {
    p <- Seurat::DimPlot(
      object,
      cols = color_list,
      label = TRUE,
      label.size = 3,
      reduction = reduction,
      group.by = group_by
    ) + Seurat::NoLegend()
  }

  if (integrate) {
    object_integrate <- Seurat::SplitObject(object, split.by = group_by)
    features <- Seurat::SelectIntegrationFeatures(object.list = object_integrate)
    object_integrate <- Seurat::FindIntegrationAnchors(object.list = object_integrate, anchor.features = features)
    object_integrate <- Seurat::IntegrateData(anchorset = object_integrate, features = features)
    object[["integrated"]] <- object_integrate@assays$integrated

    Seurat::DefaultAssay(object) <- "integrated"
    object <- Seurat::FindVariableFeatures(object)
    object <- Seurat::ScaleData(object)
    object <- Seurat::RunPCA(object)
    object <- Seurat::FindNeighbors(object, dims = dims_num)
    object <- Seurat::FindClusters(object, resolution = resolution)
    object <- Seurat::RunUMAP(object, dims = dims_num)

    if (plot) {
      p <- p + Seurat::DimPlot(
        object,
        cols = color_list,
        label = TRUE,
        label.size = 3,
        reduction = reduction,
        group.by = group_by
      ) + Seurat::NoLegend()
    }
  }

  if (plot) {
    print(p)
  }

  return(object)
}

barplot_function <- function(enrichr, font.size = NULL) {
  enrich_results <- enrichr[[1]]
  enrich_results$Term <- factor(
    enrich_results$Term,
    levels = enrich_results$Term[order(enrich_results$Combined.Score)]
  )
  enrich_results <- enrich_results[1:10, ]

  g <- ggplot(
    enrich_results,
    aes(
      x = Term,
      y = Combined.Score,
      fill = P.value, # Adjusted.P.value
      color = P.value
    )
  ) +
    geom_bar(
      stat = "identity"
    ) + coord_flip() +
    # scale_fill_distiller(palette = "Blues") +
    # scale_fill_distiller(palette = "Greys") +
    labs(y = "Combined Score", x = "") +
    theme_bw()
  g
}

enrich_functions <- function(
    genes_list,
    method = "clusterProfiler",
    enrich_method = "GO",
    plot = TRUE,
    plot_method = "dotplot",
    figure_save = FALSE,
    filename = NULL) {
  require(org.Hs.eg.db)
  require(ggplot2)

  method <- match.arg(method, c("clusterProfiler", "enrichR"))
  enrich_method <- match.arg(enrich_method, c("GO", "KEGG", "CellType"))

  # if (enrich_method == "GO" || enrich_method == "KEGG") {
  #   method <- "clusterProfiler"
  #   barplot <- graphics::barplot
  # }
  # if (enrich_method == "CellType") {
  #   method <- "enrichR"
  #
  #   enrich_res <- enrichR::enrichr(genes_list, databases = "Azimuth_Cell_Types_2021")
  #   # enrich_res <- enrichR::enrichr(genes_list, databases = "KEGG_2021_Human")
  #   # enrich_res <- enrichR::enrichr(genes_list, databases = "GO_Molecular_Function_2023")
  #   plot_method <- "barplot"
  #   barplot <- barplot_function
  # }
  if (method == "enrichR") {
    library(enrichR)
    enrich_res <- switch (
      enrich_method,
      "CellType" = enrichR::enrichr(genes_list, databases = "Azimuth_Cell_Types_2021"),
      "KEGG" = enrichR::enrichr(genes_list, databases = "KEGG_2021_Human"),
      "GO" = enrichR::enrichr(genes_list, databases = "GO_Biological_Process_2023") # Biological_Process # Cellular_Component # Molecular_Function
    )
    plot_method <- "barplot"
    barplot <- barplot_function
  }
  if (method == "clusterProfiler") {
    barplot <- graphics::barplot

    gene <- clusterProfiler::bitr(
      genes_list,
      fromType = 'SYMBOL',
      toType = 'ENTREZID',
      OrgDb = 'org.Hs.eg.db'
    )
    gene <- gene$ENTREZID

    if (enrich_method == "GO") {
      enrich_res <- clusterProfiler::enrichGO(
        gene = gene,
        OrgDb = "org.Hs.eg.db",
        ont = "ALL",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = FALSE
      )
    }
    if (enrich_method == "KEGG") {
      enrich_res <- clusterProfiler::enrichKEGG(
        gene = gene,
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05
      )
    }
  }

  if (plot) {
    plot_method <- match.arg(
      plot_method,
      c("dotplot", "barplot", "cnetplot", "emapplot", "heatplot")
    )

    p <- switch (
      plot_method,
      "dotplot" = enrichplot::dotplot(enrich_res, font.size = 15),
      "barplot" = barplot(enrich_res, font.size = 15),
      "cnetplot" = enrichplot::cnetplot(enrich_res, font.size = 15),
      "emapplot" = enrichplot::emapplot(enrichplot::pairwise_termsim(enrich_res), font.size = 15),
      "heatplot" = enrichplot::heatplot(enrich_res, showCategory = 10)
    )
    print(p)

    if (figure_save) {
      if (is.null(filename)) {
        filename <- paste0(paste(method, enrich_method, plot_method, sep = "_"), ".pdf")
      }
      ggsave(filename = filename, p, width = 8, height = 6)
    }
  }
}
