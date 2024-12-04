#' @include setClass.R
#' @include setGenerics.R

#' @title Fit models for gene expression
#'
#' @description Fit statistical models to infer gene regulatory relationships
#'
#' @param object The input data for model fitting
#' @param regulators A character vector of regulators
#' @param targets A character vector of target genes
#' @param network_name Name for the inferred network
#' @param peak_to_gene_method Method to link peaks to genes ("Signac" or "GREAT")
#' @param upstream Distance upstream of TSS to consider
#' @param downstream Distance downstream of TSS to consider
#' @param extend Extension distance for GREAT method
#' @param only_tss Whether to measure distance from TSS only
#' @param peak_to_gene_domains Custom regulatory domains
#' @param gene_cor_threshold Correlation threshold for TFs
#' @param peak_cor_threshold Correlation threshold for peaks
#' @param aggregate_rna_col Column name for RNA aggregation
#' @param aggregate_peaks_col Column name for peaks aggregation
#' @param method Statistical method for model fitting
#' @param interaction_term Type of interaction variable
#' @param adjust_method Method for p-value adjustment
#' @param scale Whether to scale the data
#' @param verbose Whether to show progress messages
#' @param cores Number of cores for parallel processing
#' @param celltype Current cell type being processed
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing model coefficients and fit statistics
#'
#' @md
#' @docType methods
#' @rdname fit_models
#' @export
setGeneric(
  name = "fit_models",
  signature = c("object"),
  def = function(object,
                 regulators = NULL,
                 targets = NULL,
                 network_name = paste0(method, "_network"),
                 peak_to_gene_method = c("Signac", "GREAT"),
                 upstream = 100000,
                 downstream = 0,
                 extend = 1000000,
                 only_tss = FALSE,
                 peak_to_gene_domains = NULL,
                 gene_cor_threshold = 0.1,
                 peak_cor_threshold = 0.1,
                 aggregate_rna_col = NULL,
                 aggregate_peaks_col = NULL,
                 method = c(
                   "srm",
                   "glm",
                   "glmnet",
                   "cv.glmnet",
                   "brms",
                   "xgb"
                 ),
                 interaction_term = ":",
                 adjust_method = "fdr",
                 scale = FALSE,
                 verbose = TRUE,
                 cores = 1,
                 celltype = NULL,
                 ...) {
    UseMethod(
      generic = "fit_models",
      object = object
    )
  }
)

#' @rdname fit_models
#' @export
setMethod(
  f = "fit_models",
  signature = signature(object = "CSNObject"),
  definition = function(object,
                        regulators = NULL,
                        targets = NULL,
                        network_name = paste0(method, "_network"),
                        peak_to_gene_method = c("Signac", "GREAT"),
                        upstream = 100000,
                        downstream = 0,
                        extend = 1000000,
                        only_tss = FALSE,
                        peak_to_gene_domains = NULL,
                        gene_cor_threshold = 0.1,
                        peak_cor_threshold = 0.1,
                        aggregate_rna_col = NULL,
                        aggregate_peaks_col = NULL,
                        method = c(
                          "srm",
                          "glm",
                          "glmnet",
                          "cv.glmnet",
                          "brms",
                          "xgb"
                        ),
                        interaction_term = ":",
                        adjust_method = "fdr",
                        scale = FALSE,
                        verbose = TRUE,
                        cores = 1,
                        celltype = NULL,
                        ...) {
    method <- match.arg(method)
    params <- Params(object)

    has_peaks <- !is.null(params$peak_assay)
    if (has_peaks) {
      peak_to_gene_method <- match.arg(peak_to_gene_method)
      motif2tf <- NetworkTFs(object)
      if (is.null(motif2tf)) {
        stop(
          "Motif matches have not been found. ",
          "Please run `find_motifs()` first."
        )
      }
      gene_annot <- Signac::Annotation(
        GetAssay(object, params$peak_assay)
      )
      if (is.null(gene_annot)) {
        stop("Please provide a gene annotation for the ChromatinAssay.")
      }
    } else {
      log_message(
        "No peak data found. Running RNA-only network inference.",
        verbose = verbose
      )
    }

    celltype_cells <- get_attribute(
      object,
      celltypes = celltype,
      attribute = "cells"
    )

    if (is.null(aggregate_rna_col)) {
      gene_data <- Matrix::t(
        LayerData(
          object,
          assay = params$rna_assay,
          layer = "data"
        )[, celltype_cells]
      )
      gene_groups <- TRUE
    } else {
      gene_data <- GetAssaySummary(
        object,
        assay = params$rna_assay,
        group_name = aggregate_rna_col,
        verbose = FALSE
      )
      gene_groups <- object@data@meta.data[[aggregate_rna_col]]
    }

    if (is.null(targets)) {
      targets <- get_attribute(
        object,
        celltypes = celltype,
        attribute = "genes"
      )
    }
    if (has_peaks) {
      celltype_peaks <- get_attribute(
        object,
        celltypes = celltype,
        attribute = "peaks"
      )

      if (is.null(aggregate_peaks_col)) {
        peak_data <- Matrix::t(
          LayerData(
            object,
            assay = params$peak_assay,
            layer = "data"
          )[, celltype_cells]
        )
        peak_groups <- TRUE
      } else {
        peak_data <- GetAssaySummary(
          object,
          assay = params$peak_assay,
          group_name = aggregate_peaks_col,
          verbose = FALSE
        )
        peak_groups <- object@data@meta.data[[aggregate_peaks_col]]
      }

      features <- intersect(gene_annot$gene_name, targets) |>
        intersect(rownames(GetAssay(object, params$rna_assay)))
      gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]

      celltype_ranges <- Signac::StringToGRanges(celltype_peaks)

      regions <- NetworkRegions(object)
      if (is.list(regions@motifs)) {
        celltype_motifs <- regions@motifs[[celltype]]$motifs
        if (!is.null(celltype_motifs)) {
          peaks2motif <- celltype_motifs@data
          peak_data <- peak_data[, celltype_peaks]
          colnames(peak_data) <- celltype_peaks
        } else {
          stop("No motifs data found for celltype ", celltype)
        }
      } else {
        peaks2motif <- regions@motifs@data
        peak_data <- peak_data[, celltype_peaks]
        colnames(peak_data) <- celltype_peaks
      }

      if (is.null(peak_to_gene_domains)) {
        log_message(
          "Selecting candidate regulatory regions near genes",
          verbose = verbose
        )
        peaks_near_gene <- find_peaks_near_genes(
          peaks = celltype_ranges,
          method = peak_to_gene_method,
          genes = gene_annot,
          upstream = upstream,
          downstream = downstream,
          only_tss = only_tss
        )
      } else {
        log_message(
          "Selecting candidate regulatory regions in provided domains",
          verbose = verbose
        )
        peaks_near_gene <- find_peaks_near_genes(
          peaks = celltype_ranges,
          method = "Signac",
          genes = peak_to_gene_domains,
          upstream = 0,
          downstream = 0,
          only_tss = FALSE
        )
      }

      peaks2gene <- aggregate_matrix(
        t(peaks_near_gene),
        groups = colnames(peaks_near_gene),
        fun = "sum"
      )

      peaks_at_gene <- as.logical(
        col_maxs(peaks2gene)
      )
      peaks_with_motif <- as.logical(
        row_maxs(peaks2motif * 1)
      )

      peaks_use <- peaks_at_gene & peaks_with_motif
      peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
      peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
      peak_data <- peak_data[, peaks_use, drop = FALSE]

      log_message(
        "Preparing model input",
        verbose = verbose
      )
      regulators <- intersect(
        regulators %s% colnames(motif2tf),
        colnames(motif2tf)
      )
      motif2tf <- motif2tf[, regulators, drop = FALSE]

      log_message(
        sprintf(
          "Fitting models for %d regulators and %d target genes",
          length(regulators), length(features)
        ),
        verbose = verbose
      )
      names(features) <- features
      model_fits <- parallelize_fun(
        features, function(g) {
          # Select peaks near gene
          if (!g %in% rownames(peaks2gene)) {
            log_message(
              g, " not found in 'EnsDb'",
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }
          gene_peaks <- as.logical(peaks2gene[g, ])
          if (sum(gene_peaks) == 0) {
            log_message(
              "no peaks found near ", g,
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }

          g_x <- gene_data[gene_groups, g, drop = FALSE]
          peak_x <- peak_data[peak_groups, gene_peaks, drop = FALSE]

          if (sum(colSums(peak_x != 0) >= 3) == 0) {
            log_message(
              sprintf("not enough non-zero values in peaks for %s", g),
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }

          peak_g_cor <- as(
            sparse_cor(peak_x, g_x),
            "generalMatrix"
          )
          peak_g_cor[is.na(peak_g_cor)] <- 0
          peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > peak_cor_threshold]

          if (length(peaks_use) == 0) {
            log_message(
              sprintf("no correlating peaks found for %s", g),
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }

          peak_x <- peak_x[, peaks_use, drop = FALSE]
          peak_motifs <- peaks2motif[peaks_use, , drop = FALSE]

          gene_peak_tfs <- purrr::map(
            rownames(peak_motifs), function(p) {
              x <- as.logical(peak_motifs[p, ])
              peak_tfs <- col_maxs(
                motif2tf[x, , drop = FALSE]
              )
              peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
              peak_tfs <- intersect(peak_tfs, colnames(gene_data))
              peak_tfs <- setdiff(peak_tfs, g)
              return(peak_tfs)
            }
          )
          names(gene_peak_tfs) <- rownames(peak_motifs)

          # Check correlation of peaks with target gene
          gene_tfs <- purrr::reduce(gene_peak_tfs, union)
          tf_x <- gene_data[gene_groups, gene_tfs, drop = FALSE]
          tf_g_cor <- as(sparse_cor(tf_x, g_x), "generalMatrix")
          tf_g_cor[is.na(tf_g_cor)] <- 0
          tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > gene_cor_threshold]
          if (length(tfs_use) == 0) {
            log_message(
              sprintf("no correlating TFs found for %s", g),
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }
          tf_g_corr_df <- tibble::as_tibble(
            tf_g_cor[unique(tfs_use), , drop = FALSE],
            rownames = "tf",
            .name_repair = "check_unique"
          ) |>
            dplyr::rename("tf" = 1, "corr" = 2)

          # Filter TFs and make formula string
          frml_string <- purrr::map(names(gene_peak_tfs), function(p) {
            peak_tfs <- gene_peak_tfs[[p]]
            peak_tfs <- peak_tfs[peak_tfs %in% tfs_use]
            if (length(peak_tfs) == 0) {
              return()
            }
            peak_name <- stringr::str_replace_all(p, "-", "_")
            tf_name <- stringr::str_replace_all(peak_tfs, "-", "_")
            formula_str <- paste(
              paste(peak_name, interaction_term, tf_name, sep = " "),
              collapse = " + "
            )
            return(list(tfs = peak_tfs, frml = formula_str))
          })
          frml_string <- frml_string[!purrr::map_lgl(frml_string, is.null)]
          if (length(frml_string) == 0) {
            log_message(
              sprintf("no valid peak:TF pairs found for %s", g),
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          }

          target <- stringr::str_replace_all(g, "-", "_")
          model_frml <- stats::as.formula(
            paste0(
              target, " ~ ",
              paste0(
                purrr::map(
                  frml_string,
                  function(x) x$frml
                ),
                collapse = " + "
              )
            )
          )

          nfeats <- sum(
            purrr::map_dbl(
              frml_string,
              function(x) length(x$tfs)
            )
          )
          gene_tfs <- purrr::reduce(
            purrr::map(frml_string, function(x) x$tfs),
            union
          )
          gene_x <- gene_data[gene_groups, union(g, gene_tfs), drop = FALSE]
          model_mat <- as.data.frame(cbind(gene_x, peak_x))
          if (scale) {
            model_mat <- as.data.frame(scale(as_matrix(model_mat)))
          }
          colnames(model_mat) <- stringr::str_replace_all(
            colnames(model_mat), "-", "_"
          )

          log_message(
            "Fitting model with ", nfeats, " variables for ", g,
            verbose = verbose == 2
          )
          result <- try(
            fit_model(
              model_frml,
              data = model_mat,
              method = method,
              ...
            ),
            silent = TRUE
          )
          if (any(class(result) == "try-error")) {
            log_message(
              sprintf("fitting model failed for %s", g),
              verbose = verbose == 2,
              message_type = "warning"
            )
            log_message(
              result,
              verbose = verbose == 2,
              message_type = "warning"
            )
            return()
          } else {
            result$metrics$nvariables <- nfeats
            result$corr <- tf_g_corr_df
            return(result)
          }
        },
        verbose = verbose, cores = cores
      )
    } else {
      if (is.null(regulators)) {
        regulators <- get_attribute(
          object,
          celltypes = celltype,
          attribute = "genes"
        )
      }
      model_fits <- fit_models(
        object = as_matrix(gene_data),
        regulators = regulators,
        targets = targets,
        gene_cor_threshold = gene_cor_threshold,
        method = method,
        scale = scale,
        cores = cores,
        verbose = verbose,
        ...
      )
    }

    model_fits <- model_fits[!purrr::map_lgl(model_fits, is.null)]
    if (length(model_fits) == 0) {
      log_message(
        "fitting model failed for all genes.",
        verbose = verbose == 2,
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
      variable = if (has_peaks) interaction_term else NULL,
      adjust_method = adjust_method
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
    metrics <- purrr::map_dfr(
      model_fits,
      function(x) x$metrics,
      .id = "target"
    )

    params <- list()
    params[["method"]] <- method
    params[["family"]] <- family
    params[["gene_cor_threshold"]] <- gene_cor_threshold
    if (has_peaks) {
      params[["dist"]] <- c(
        "upstream" = upstream,
        "downstream" = downstream
      )
      params[["only_tss"]] <- only_tss
      params[["interaction"]] <- interaction_term
      params[["peak_cor_threshold"]] <- peak_cor_threshold
    }

    network_obj <- methods::new(
      Class = "Network",
      regulators = regulators,
      targets = targets,
      coefficients = coefficients,
      fit = metrics,
      params = params
    )
    object@networks[[network_name]][[celltype]] <- network_obj
    object@active_network <- network_name

    return(object)
  }
)

#' @rdname fit_models
#' @export
setMethod(
  f = "fit_models",
  signature = signature(object = "matrix"),
  definition = function(object,
                        regulators = NULL,
                        targets = NULL,
                        gene_cor_threshold = 0.1,
                        method = c(
                          "srm",
                          "glm",
                          "glmnet",
                          "cv.glmnet",
                          "brms",
                          "xgb"
                        ),
                        scale = FALSE,
                        verbose = TRUE,
                        cores = 1,
                        ...) {
    method <- match.arg(method)
    gene_data <- object

    names(targets) <- targets
    log_message(
      sprintf(
        "Fitting models for %d regulators and %d target genes",
        length(regulators), length(targets)
      ),
      verbose = verbose
    )

    model_fits <- parallelize_fun(
      targets, function(g) {
        regulators_use <- setdiff(regulators, g)
        g_x <- gene_data[, g, drop = FALSE]
        tf_x <- gene_data[, regulators_use, drop = FALSE]

        tf_g_cor <- as(sparse_cor(tf_x, g_x), "generalMatrix")
        tf_g_cor[is.na(tf_g_cor)] <- 0
        tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > gene_cor_threshold]

        if (length(tfs_use) == 0) {
          log_message(
            sprintf("no correlating TFs found for %s", g),
            verbose = verbose == 2,
            message_type = "warning"
          )
          return(NULL)
        }

        tf_g_corr_df <- tibble::as_tibble(
          tf_g_cor[tfs_use, , drop = FALSE],
          rownames = "tf",
          .name_repair = "check_unique"
        ) |>
          dplyr::rename("tf" = 1, "corr" = 2)

        model_mat <- as.data.frame(tf_x[, tfs_use, drop = FALSE])
        model_mat[[g]] <- as.vector(g_x)

        if (scale) {
          model_mat <- as.data.frame(scale(model_mat))
        }

        colnames(model_mat) <- stringr::str_replace_all(
          colnames(model_mat), "-", "_"
        )

        target <- stringr::str_replace_all(g, "-", "_")
        model_frml <- stats::as.formula(
          paste0(
            target, " ~ ",
            paste(
              stringr::str_replace_all(tfs_use, "-", "_"),
              collapse = " + "
            )
          )
        )

        result <- try(
          fit_model(
            model_frml,
            data = model_mat,
            method = method,
            ...
          ),
          silent = TRUE
        )

        if (inherits(result, "try-error")) {
          log_message(
            sprintf("fitting model failed for %s", g),
            verbose = verbose == 2,
            message_type = "warning"
          )
          return(NULL)
        }

        result$corr <- tf_g_corr_df
        return(result)
      },
      verbose = verbose,
      cores = cores
    )

    return(model_fits)
  }
)
