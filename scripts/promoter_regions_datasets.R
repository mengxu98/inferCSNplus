prepare.data <- function(hg) {
  promoter_regions <- read.table(file = paste0("../inferCSN_data/", hg, ".promoter.regions.txt"))
  names(promoter_regions) <- c("Chrom","Starts","Ends","genes")
  genes <- lapply(promoter_regions$genes, function(x) strsplit(x,"[|]")[[1]][1])
  genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
  genes <- unlist(genes)
  promoter_regions$genes <- genes
  unik <- !duplicated(genes)
  promoter_regions <- promoter_regions[unik,]
  return(promoter_regions)
}

promoter_regions_hg19 <- prepare.data("hg19")
promoter_regions_hg38 <- prepare.data("hg38")

usethis::use_data(promoter_regions_hg19)
usethis::use_data(promoter_regions_hg38)
