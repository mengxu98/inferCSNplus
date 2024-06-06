source("scripts/functions.R")
library(Seurat)

Ma_Sestan <- readRDS("~/HAR/brain_data/Ma_Sestan_mat.rds")
seurat_object <- seurat_function(
  Ma_Sestan,
  re_create_object = TRUE
)
