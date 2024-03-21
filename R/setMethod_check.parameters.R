#' @title Check input parameters
#'
#' @param matrix An expression matrix, cells by genes
#' @inheritParams inferCSN
#'
#' @return No return value, called for check input parameters
#' @export
check.parameters <- function(
    matrix,
    penalty,
    algorithm,
    cross_validation,
    seed,
    n_folds,
    k_folds,
    r_threshold,
    regulators,
    targets,
    regulators_num,
    verbose,
    cores) {
  if (verbose) message("Checking input parameters.")

  error_message <- paste0(
    "Parameter matrix must be a two-dimensional matrix,
    where each column corresponds to a gene and each row corresponds to a sample/cell."
  )
  if (!is.matrix(matrix) && !is.array(matrix)) {
    stop(error_message)
  }

  if (length(dim(matrix)) != 2) {
    stop(error_message)
  }

  if (is.null(colnames(matrix))) {
    stop("Parameter matrix must contain the names of the genes as colnames.")
  }

  # Check the penalty term of the regression model
  if (!any(c("L0", "L0L2") == penalty)) {
    stop(
      "inferCSN does not support ", penalty, " penalty regression.\n",
      "Please set penalty item as 'L0' or 'L0L2'."
    )
  }

  # Check the algorithm of the regression model
  if (!any(c("CD", "CDPSI") == algorithm)) {
    stop(
      "inferCSN does not support ", algorithm, " algorithmn.\n",
      "Please set algorithm item as 'CD' or 'CDPSI'."
    )
  }

  if (!is.numeric(seed)) {
    seed <- 1
    warning("Supplied seed is not a valid integer, initialize 'seed' to 1.")
  }

  if (!is.null(k_folds)) {
    if (!(is.numeric(k_folds) && k_folds > 0 && k_folds < 10)) {
      stop("Please set 'k_folds' value between: (0, 10].")
    }
  }

  if (!is.null(targets)) {
    if (!is.vector(targets)) {
      stop("Parameter 'targets' must a vector (of indices or gene names).")
    }

    if (is.numeric(targets)) {
      if (max(targets) > nrow(matrix)) stop("At least one index in 'targets' exceeds the number of genes.")
      if (min(targets) <= 0) stop("The indexes in 'targets' should be >=1.")
    }

    if (any(table(targets) > 1)) stop("Please provide each target only once.")

    if (is.character(targets)) {
      targetsInMatrix <- intersect(targets, colnames(matrix))
      if (length(targetsInMatrix) == 0) {
        stop("The genes must contain at least one target.")
      }

      if (length(targetsInMatrix) < length(targets)) {
        warning("Only ", length(targetsInMatrix), " out of ", length(targets), " candidate regulators are in the expression matrix.")
      }
    }
  }

  if (!is.null(regulators)) {
    if (is.list(regulators)) {
      if (!all(names(regulators) %in% targets)) {
        stop("Regulators: If provided as a named list, all names should be targets.")
      }
      regulators <- unique(unlist(regulators))
    }
    if (!is.null(regulators)) {
      if (length(regulators) < 2) stop("Provide at least 2 potential regulators.")

      if (!is.vector(regulators)) {
        stop("Parameter 'regulators' must a vector of indices or gene names.")
      }

      if (is.numeric(regulators)) {
        if (max(regulators) > nrow(matrix)) {
          stop("At least one index in 'regulators' exceeds the number of genes.")
        }
        if (min(regulators) <= 0) stop("The indexes in 'regulators' should be >=1.")
      }

      if (any(table(regulators) > 1)) stop("Please provide each regulator only once.")

      if (is.character(regulators)) {
        regulatorsInMatrix <- intersect(regulators, colnames(matrix))
        if (length(regulatorsInMatrix) < 2) {
          stop("Fewer than 2 regulators in the columns of expression matrix.")
        }

        if (length(regulatorsInMatrix) < length(regulators)) {
          warning("Only ", length(regulatorsInMatrix), " out of ", length(regulators), " candidate regulators are in the expression matrix.")
        }
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    stop(x = "Parameter cores should be a stricly positive integer.")
  }

  if (verbose) message("All parameters check done.")

  if (verbose) message("Using '", penalty, "' penalty.")
  if (verbose && cross_validation) message("Using cross validation.")
}
