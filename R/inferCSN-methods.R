#' @include setClass.R
#' @include setGenerics.R
#' @include inferCSN.R

#' @rdname inferCSN
#' @method inferCSN Network
#' @export
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' object <- initiate_object(example_matrix)
#' object <- inferCSN(object)
#' }
setMethod(
  f = "inferCSN",
  signature = signature(object = "Network"),
  definition = function(object,
                        penalty = "L0",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 5,
                        subsampling_method = "sample",
                        subsampling_ratio = 1,
                        r_squared_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        cores = 1,
                        verbose = TRUE,
                        method = c(
                          "srm",
                          "glm",
                          "glmnet",
                          "cv.glmnet",
                          "xgb",
                          "susie"
                        ),
                        gene_cor_threshold = 0,
                        ...) {
    matrix <- thisutils::as_matrix(object@data)
    if (is.null(regulators)) {
      regulators <- object@regulators
      if (length(regulators) == 0) {
        regulators <- colnames(matrix)
      }
    } else {
      regulators <- intersect(
        regulators,
        colnames(matrix)
      )
    }
    if (is.null(targets)) {
      targets <- object@targets
      if (length(targets) == 0) {
        targets <- colnames(matrix)
      }
    } else {
      targets <- intersect(
        targets,
        colnames(matrix)
      )
    }

    model_fits <- fit_models(
      object = matrix,
      regulators = regulators,
      targets = targets,
      gene_cor_threshold = gene_cor_threshold,
      method = method,
      scale = FALSE,
      cores = cores,
      verbose = verbose,
      ...
    )
    model_fits <- model_fits[!purrr::map_lgl(model_fits, is.null)]
    if (length(model_fits) == 0) {
      thisutils::log_message(
        "fitting model failed for all genes.",
        verbose = verbose,
        message_type = "warning"
      )
    }

    coefficients <- purrr::map_dfr(
      model_fits,
      function(x) x$coefficients,
      .id = "target"
    )
    coefficients <- format_coefs(
      coefficients,
      variable = NULL,
      adjust_method = "fdr"
    )
    corrs <- purrr::map_dfr(
      model_fits,
      function(x) x$corr,
      .id = "target"
    )
    if (nrow(coefficients) > 0) {
      coefficients <- suppressMessages(
        dplyr::left_join(coefficients, corrs)
      )
    }
    object@metrics <- purrr::map_dfr(
      model_fits,
      function(x) x$metrics,
      .id = "target"
    )

    object@coefficients <- coefficients
    object@regulators <- regulators
    object@targets <- targets
    object@params <- list(
      method = method,
      penalty = penalty,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      gene_cor_threshold = gene_cor_threshold,
      r_squared_threshold = r_squared_threshold,
      cores = cores,
      verbose = verbose
    )
    object <- .process_Network(
      object,
      r_squared_threshold = r_squared_threshold
    )
    object@network <- export_csn(object)

    return(object)
  }
)

#' @param celltypes Character vector of cell types to infer networks for.
#' @param network_name network_name.
#' @param peak_to_gene_method Character specifying the method to
#' link peak overlapping motif regions to nearby genes. One of \code{Signac} or \code{GREAT}.
#' @param upstream Integer defining the distance upstream of the gene to consider as potential regulatory region.
#' @param downstream Integer defining the distance downstream of the gene to consider as potential regulatory region.
#' @param extend Integer defining the distance from the upstream and downstream of the basal regulatory region.
#' Only used of `peak_to_gene_method = 'GREAT'`.
#' @param only_tss Logical. Measure distance from the TSS (\code{TRUE}) or from the entire gene body (\code{FALSE}).
#' @param peak_to_gene_domains \code{GenomicRanges} object with regulatory regions for each gene.
#' @param gene_cor_threshold Threshold for TF - target gene correlation.
#' @param peak_cor_threshold Threshold for binding peak - target gene correlation.
#' @param aggregate_rna_col aggregate_rna_col
#' @param aggregate_peaks_col aggregate_peaks_col
#' @param method A character string indicating the method to fit the model.
#' * \code{'srm'} - Sparse Regression Model.
#' * \code{'glm'} - Generalized Liner Model with \code{\link[stats]{glm}}.
#' * \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[stats]{family}} for mode details.
#' @param interaction_term The interaction variable to use in the model between TF and binding site.
#' * \code{'+'} for additive interaction.
#' * \code{':'} for 'multiplicative' interaction.
#' * \code{'*'} for crossing interaction, i.e. additive AND 'multiplicative'.
#' For more info, see \code{\link[stats]{formula}}
#' @param scale Logical. Whether to z-transform the expression and accessibility matrices.
#' @param adjust_method Method for adjusting p-values.
#'
#' @return A CSNObject.
#'
#' @rdname inferCSN
#' @method inferCSN CSNObject
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(object = "CSNObject"),
  definition = function(object,
                        penalty = "L0",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 5,
                        subsampling_method = "sample",
                        subsampling_ratio = 1,
                        r_squared_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        cores = 1,
                        verbose = TRUE,
                        celltypes = NULL,
                        network_name = paste0(method, "_network"),
                        peak_to_gene_method = c("Signac", "GREAT"),
                        upstream = 100000,
                        downstream = 0,
                        extend = 1000000,
                        only_tss = FALSE,
                        peak_to_gene_domains = NULL,
                        gene_cor_threshold = 0.1,
                        peak_cor_threshold = 0.,
                        aggregate_rna_col = NULL,
                        aggregate_peaks_col = NULL,
                        method = c(
                          "srm",
                          "glm",
                          "glmnet",
                          "cv.glmnet",
                          "xgb",
                          "susie"
                        ),
                        alpha = 0.5,
                        family = "gaussian",
                        interaction_term = ":",
                        adjust_method = "fdr",
                        scale = FALSE,
                        ...) {
    method <- match.arg(method)
    peak_to_gene_method <- match.arg(peak_to_gene_method)

    if (is.null(celltypes)) {
      celltypes <- object@metadata$celltypes
    }

    for (celltype in celltypes) {
      thisutils::log_message(
        "Inferring network for celltype: '", celltype, "'",
        verbose = verbose
      )
      celltype_genes <- get_attribute(
        object,
        celltypes = celltype,
        attribute = "genes"
      )

      celltype_genes <- targets %ss% celltype_genes

      if (length(celltype_genes) == 0) {
        thisutils::log_message(
          "no target genes found for ", celltype, ", skipping",
          verbose = verbose,
          message_type = "warning"
        )
        next
      }

      object <- fit_models(
        object = object,
        regulators = regulators,
        targets = celltype_genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        gene_cor_threshold = gene_cor_threshold,
        peak_cor_threshold = peak_cor_threshold,
        aggregate_rna_col = aggregate_rna_col,
        aggregate_peaks_col = aggregate_peaks_col,
        method = method,
        alpha = alpha,
        family = family,
        interaction_term = interaction_term,
        adjust_method = adjust_method,
        scale = scale,
        verbose = verbose,
        cores = cores,
        celltype = celltype,
        ...
      )
    }

    thisutils::log_message(
      "Network inference summary:",
      verbose = verbose
    )

    object <- .process_csn(
      object,
      r_squared_threshold = r_squared_threshold
    )
    for (celltype in celltypes) {
      network <- export_csn(
        object,
        celltypes = celltype
      )
      if (!is.null(network)) {
        thisutils::log_message(
          sprintf(
            "   %s: %d edges (%d regulators -> %d targets)",
            celltype,
            nrow(network),
            length(unique(network$regulator)),
            length(unique(network$target))
          ),
          verbose = verbose,
          cli_model = FALSE
        )
      } else {
        thisutils::log_message(
          celltype, ": no edges inferred",
          verbose = verbose,
          message_type = "warning"
        )
      }
    }

    return(object)
  }
)

#' @rdname inferCSN
#' @method inferCSN Seurat
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(object = "Seurat"),
  definition = function(object, ...) {
    stop(
      "inferCSN is not supported for Seurat object",
      "\n  Please run `initiate_object()` for Seurat object first"
    )
  }
)
