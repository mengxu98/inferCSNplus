#' @include setClass.R
#' @include setGenerics.R

#' @title Calculate gene ranks in network
#' @description Calculate gene importance using PageRank algorithm
#'
#' @param object A Network object or data frame containing network information
#' @param regulators Optional character vector of regulator genes
#' @param targets Optional character vector of target genes
#' @param directed Logical, whether the network is directed
#'
#' @return A data frame of gene ranks
#' @export
setGeneric(
  name = "calculate_gene_rank",
  signature = c("object"),
  def = function(object,
                 regulators = NULL,
                 targets = NULL,
                 directed = FALSE) {
    standardGeneric("calculate_gene_rank")
  }
)

#' @rdname calculate_gene_rank
#' @export
setMethod(
  "calculate_gene_rank",
  signature(object = "Network"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- as.data.frame(object@network)
    .calculate_page_rank(network_table, directed)
  }
)

#' @rdname calculate_gene_rank
#' @export
setMethod(
  "calculate_gene_rank",
  signature(object = "data.frame"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- network_format(
      object,
      regulators,
      targets,
      abs_weight = FALSE
    )
    .calculate_page_rank(
      network_table,
      directed
    )
  }
)

.calculate_page_rank <- function(network_table, directed) {
  network <- igraph::graph_from_data_frame(
    network_table,
    directed = directed
  )
  page_rank_res <- data.frame(
    igraph::page_rank(network, directed = directed)$vector
  )
  colnames(page_rank_res) <- c("rank_weight")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "rank_weight")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$rank_weight,
    decreasing = TRUE
  ), ]
  page_rank_res$regulator <- ifelse(
    page_rank_res$gene %in% unique(network_table$regulator),
    "TRUE", "FALSE"
  )
  rownames(page_rank_res) <- NULL

  return(page_rank_res)
}
