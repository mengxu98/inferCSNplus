#' run inferCSN on Seurat object
#'
#' @param object Seurat object.
#' @param peakcalling  call peak
#' @param macs2.path  path to macs2
#' @param fragments  fragments file
#' @param k_neigh Number of cells to aggregate per group.
#' @param atacbinary Logical, should accessibility values be binarized
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param reduction.name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, whether need to do size normalization
#' @param genome.info the TSS information of genome, e.g. hg19, hg38
#' @param focus_markers the focused genes
#' @param params the list of parameters used in Xgboost
#' @param cores  the number of threads can be manually specified in Xgboost trainning stage, default is 2
#' @param early_stop Logical, whether use early stop rule on validation data to reduce overfitting
#' @param HC_cutoff the threshold of high functional CREs
#' @param LC_cutoff the threshold of low functional CREs
#' @param rescued Logical, whether to rescue highly correlated CREs
#' @param seed Random seed
#' @param verbose Logical, should warning and info messages be printed
#' @import Signac
#' @import Seurat
#' @import Matrix
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom cicero find_overlapping_coordinates
#' @return a Seurat object with new links assay.
#' @export
data.processing <- function(object,
                            aggregate = TRUE,
                            peakcalling = FALSE,
                            macs2.path = NULL,
                            fragments = NULL,
                            k_neigh = 50,
                            atacbinary = TRUE,
                            max_overlap = 0.8,
                            reduction.name = NULL,
                            size_factor_normalize = FALSE,
                            genome.info,
                            focus_markers = NULL,
                            cores = 2,
                            early_stop = FALSE, # crossvaliation + sample?
                            HC_cutoff = NULL,
                            LC_cutoff = NULL,
                            rescued = FALSE,
                            seed = 123,
                            verbose = TRUE) {


  options(stringsAsFactors = FALSE)
  

  inferCSN_results <- list()
  TXs <- list()
  TYs <- list()
  for (i in 1:length(focus_markers)) {
    if (verbose) {
      message("Inferring links for: ", focus_markers[i])
    }

    p1 <- paste(chr[i], ":", Starts[i] - 500, "-", Starts[i], sep = "")
    p2 <- paste(chr[i], ":", Starts[i] - 250000, "-", Starts[i] + 250000, sep = "")
    promoters <- cicero::find_overlapping_coordinates(peaks, p1)
    enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
    enhancers <- setdiff(enhancers, promoters)

    ###build data matrix of each gene used in model
    if ("RNA" %in% names(agg_data)) {
      idx <- which(rownames(data_rna) == focus_markers[i])
    } else {
      idx <- 1
    }

    if ((length(promoters) > 0 && length(enhancers) > 1) && length(idx) != 0) {
      id1 <- match(promoters, peaks)
      id1 <- id1[!is.na(id1)]
      id2 <- match(enhancers, peaks)
      id2 <- id2[!is.na(id2)]
      id2_new <- setdiff(id2, id1)
      X <- data_atac[id2_new, ]
      TXs[[i]] <- X
      Y <- data_atac[id1, ]
      if (length(id1) > 1) {
        Y <- colSums(Y)
      }
      Y <- t(as.matrix(Y))
      rownames(Y) <- peaks[id1[1]]
      TYs[[i]] <- Y
      if ("RNA" %in% names(agg_data)) {
        Z <- data_rna[idx, ]
        Z <- t(as.matrix(Z))
        rownames(Z) <- focus_markers[i]
      } else {
        Z <- Y
      }
      flag <- 1
    } else {
      flag <- 0
      message("There are less than two peaks detected within 500 kb for ", focus_markers[i])
    }

  }

  return(data_object)
}
