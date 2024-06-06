# Human microglial state dynamics in Alzheimer’s disease progression
# https://doi.org/10.1016/j.cell.2023.08.037
# data: https://compbio.mit.edu/microglia_states/
# https://cells.ucsc.edu/?ds=rosmap-ad-aging-brain+microglia-states
# https://www.dropbox.com/scl/fo/100k2i24hbapxjk34k69e/AOJwl1MIJSVh2XIp-QS806E?rlkey=k1vhpx1r7175up5fsyyavdlcb&e=1&dl=0
# code: https://github.com/nasunmit/scMicroglia

rna <- readRDS("~/HAR/brain_data/GSE227222/GSE227222_iMGs.snRNAseq.Jan2023.counts.rds")
rna_metadata <- readRDS("~/HAR/brain_data/GSE227222/GSE227222_iMGs.snRNAseq.Jan2023.metadata.rds")
rna2 <- readRDS("~/HAR/brain_data/GSE227222/GSE227222_iMGs.snRNAseq.May2022.counts.rds")
rna_metadata2 <- readRDS("~/HAR/brain_data/GSE227222/GSE227222_iMGs.snRNAseq.May2022.metadata.rds")
