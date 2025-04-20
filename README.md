# ***inferCSN+*** <img src="man/figures/logo.svg" align="right" width="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/inferCSN)](https://github.com/cran/inferCSN) [![develop-ver](https://img.shields.io/github/r-package/v/mengxu98/inferCSNplus?label=develop-ver)](https://github.com/mengxu98/inferCSNplus/blob/main/DESCRIPTION) [![R-CMD-check](https://github.com/mengxu98/inferCSNplus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mengxu98/inferCSN/actions/workflows/R-CMD-check.yaml) [![test-coverage](https://github.com/mengxu98/inferCSNplus/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/mengxu98/inferCSNplus/actions/workflows/test-coverage.yaml) [![pkgdown](https://github.com/mengxu98/inferCSNplus/actions/workflows/pkgdown.yaml/badge.svg)](https://mengxu98.github.io/inferCSNplus/reference/index.html) [![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/inferCSN)](https://CRAN.R-project.org/package=inferCSN) [![code-size](https://img.shields.io/github/languages/code-size/mengxu98/inferCSNplus)](https://github.com/mengxu98/inferCSNplus)

<!-- badges: end -->

## **Introduction**

[*`inferCSN+`*](https://mengxu98.github.io/inferCSNplus/) is an R package for ***infer***ring ***C***ell-***S***pecific gene regulatory ***N***etwork from single-cell omics data.

<img src="man/figures/inferCSNplus.svg" width="90%"/>

## **Installation**

You can install the released version from [*`CRAN`*](https://github.com/cran) use:

``` r
install.packages("inferCSN")
# or
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("inferCSN")
```

You can install the development version from [*`GitHub`*](https://github.com/mengxu98/inferCSNplus) use [*`pak`*](https://github.com/r-lib/pak):

``` r
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/inferCSNplus")
```

## **Usage**

### **Examples**

#### For a `matrix` object.

``` r
library(inferCSN)
data("example_matrix")

network <- inferCSN(
  example_matrix
)
```

#### For a `matrix` object, initiate it to `Network` object.

``` r
library(inferCSN)
data("example_matrix")

object <- initiate_object(
  example_matrix
)
object
object <- inferCSN(
  object
)

network_table <- export_csn(
  object
)
head(network_table)
```

#### For a `Seurat` object, it contains multiome data (scRNA-seq and/or scATAC-seq)

Multiome data from [GSE204684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204684).
``` r
library(inferCSN)
library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)

# initiate Seurat object to CSNObject
object <- initiate_object(
  object,
  peak_assay = "ATAC",
  filter_by = "celltype",
  celltype = NULL,
  celltype_by = "celltype"
)
object <- find_motifs(
  object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
object <- inferCSN(
  object,
  cores = 4
)
networks <- export_csn(
  object
)
head(networks[[1]])

# network processing
object <- find_modules(
  object,
  p_thresh = 0.3,
  nvar_thresh = 10,
  min_genes_per_module = 5,
  rsq_thresh = 0.05
)

modules <- NetworkModules(object)
modules[[1]]@meta

# network visualization
p <- plot_gof(object)
p[[1]]

p2 <- plot_module_metrics(object)
p2[[1]]

object <- get_network_graph(object)
plot_network_graph(object)

targets <- get_attribute(
  object,
  attribute = "targets"
)
targets <- targets[["RG"]]
tf <- targets[1:100]
object <- get_tf_network(
  object,
  tfs = tf,
  celltypes = "RG",
  verbose = FALSE
)

plot_list <- plot_tf_network(
  object,
  tfs = tf,
  edge_width = 0.5,
  celltypes = "RG",
  verbose = FALSE
)

plot_list$RG
```

More functions and usages about [*`inferCSN+`*](https://mengxu98.github.io/inferCSNplus/)? Please reference [*`here`*](https://mengxu98.github.io/inferCSNplus/reference/index.html).

## **Citation**

If you use [*`inferCSN+`*](https://github.com/mengxu98plus/inferCSN) in your work, please cite it reference [*`here`*](https://mengxu98.github.io/inferCSNplus/authors.html#citation).
