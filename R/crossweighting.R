#' Perform crossweighting
#'
#' @param weight_table GRN dataframe, the result of running reconstructargetRN or reconstructargetRN_GENIE3
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object result of running findDynGenes
#' @param lag lag window on which to run cross-correlation. Cross-correlaiton computed from -lag to +lag.
#' @param min minimum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) less than min will not be negatively weighted.
#' @param max maximum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) greater than max will have weights set to 0.
#' @param symmetric_filter whether or not to employ a symmetric weight scheme. If true, absolute offset is used in place of offset.
#' @param filter_thresh after crossweighting, edges with weights less than filter_thresh will be set to 0.
#'
#' @return weight_table with offset and weighted_score added
#'
#' @export
crossweight <- function(
    weight_table,
    matrix,
    dynamic_object,
    lag = floor(ncol(matrix) / 5),
    min = ceiling(ncol(matrix) / 50),
    max = floor(ncol(matrix) / 12),
    symmetric_filter = FALSE,
    filter_thresh = 0) {
  # order matrix
  matrix <- matrix[, rownames(dynamic_object$cells)]

  weight_table$target <- as.character(weight_table$target)
  weight_table$regulator <- as.character(weight_table$regulator)

  offset <- apply(weight_table, 1, cross_corr, matrix = matrix, lag = lag)
  weight_table$offset <- offset

  weighted_score <- c()
  for (i in 1:nrow(weight_table)) {
    new <- score_offset(
      # weight_table$zscore[i],
      weight_table$weight[i],
      weight_table$offset[i],
      min = min,
      max = max,
      symmetric_filter = symmetric_filter
    )
    weighted_score <- c(weighted_score, new)
  }

  weight_table$weighted_score <- weighted_score

  weight_table <- weight_table[weight_table$weighted_score > filter_thresh, ]

  return(weight_table)
}


cross_corr <- function(
    grn_row,
    matrix,
    lag) {
  target <- grn_row[1]
  regulator <- grn_row[2]

  x <- stats::ccf(as.numeric(matrix[regulator, ]), as.numeric(matrix[target, ]), lag, plot = FALSE)

  df <- data.frame(lag = x$lag, cor = abs(x$acf))
  df <- df[order(df$cor, decreasing = TRUE), ]
  offset <- mean(df$lag[1:ceiling((2 / 3) * lag)])

  return(offset)
}

score_offset <- function(
    score,
    offset,
    min = 2,
    max = 20,
    symmetric_filter = FALSE) {
  if (symmetric_filter) {
    offset <- abs(offset)
  }

  if (offset <= min) {
    res <- score
  } else if (offset >= max) {
    res <- 0
  } else {
    # linear weighting scheme according to y=(-x/(max-min))+1
    weight <- (-offset / (max - min)) + 1
    res <- score * weight
  }

  return(res)
}

# estimates min and max values for crossweighting
# for now assumes uniform cell density across pseudotime/only considers early time
# this needs to be refined if it's to be useful...
crossweight_params <- function(
    matrix,
    dynamic_object,
    pseudotime_min = 0.005,
    pseudotime_max = 0.01) {
  matrix <- matrix[, rownames(dynamic_object$cells)]
  ncells <- nrow(dynamic_object$cells)
  min <- nrow(dynamic_object$cells[dynamic_object$cells$pseudotime < pseudotime_min, ])
  max <- nrow(dynamic_object$cells[dynamic_object$cells$pseudotime < pseudotime_max, ])

  params <- list(min = min, max = max)

  return(params)
}
