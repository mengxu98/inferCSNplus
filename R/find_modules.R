#' @title Find TF modules in regulatory network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @rdname find_modules
#' @export find_modules
setGeneric(
  "find_modules",
  signature = "object",
  function(object, ...) {
    standardGeneric("find_modules")
  }
)

#' @param p_thresh Float indicating the significance threshold on the adjusted p-value.
#' @param rsq_thresh Float indicating the \eqn{R^2} threshold on the adjusted p-value.
#' @param nvar_thresh Integer indicating the minimum number of variables in the model.
#' @param min_genes_per_module Integer indicating the minimum number of genes in a module.
#' @param xgb_method Method to get modules from xgb models
#' \code{tf} - Choose top targets for each TF.
#' \code{target} - Choose top TFs for each target gene.
#' @param xgb_top Interger indicating how many top targets/TFs to return.
#' @param verbose Print messages.
#'
#' @return A Network object.
#'
#' @rdname find_modules
#' @export
#' @method find_modules Network
setMethod(
  f = "find_modules",
  signature = "Network",
  definition = function(object,
                        p_thresh = 0.05,
                        rsq_thresh = 0.1,
                        nvar_thresh = 10,
                        min_genes_per_module = 5,
                        xgb_method = c("tf", "target"),
                        xgb_top = 50,
                        verbose = TRUE,
                        ...) {
    fit_method <- NetworkParams(object)$method
    xgb_method <- match.arg(xgb_method)

    if (!fit_method %in% c("srm", "glm", "cv.glmnet", "glmnet", "xgb")) {
      stop(
        paste0(
          'find_modules() is not yet implemented for "',
          fit_method,
          '" models'
        )
      )
    }

    models_use <- metrics(object) |>
      dplyr::filter(r_squared > rsq_thresh & nvariables > nvar_thresh) |>
      dplyr::pull(target) |>
      unique()

    modules <- coef(object) |>
      dplyr::filter(target %in% models_use)

    if (fit_method %in% c("srm", "cv.glmnet", "glmnet")) {
      modules <- modules |>
        dplyr::filter(coefficient != 0)
    } else if (fit_method == "xgb") {
      modules <- modules |>
        dplyr::group_by_at(xgb_method) |>
        dplyr::top_n(xgb_top, gain) |>
        dplyr::mutate(coefficient = sign(corr) * gain)
    } else {
      modules <- modules |>
        dplyr::filter(ifelse(is.na(padj), T, padj < p_thresh))
    }

    modules <- modules |>
      dplyr::group_by(target) |>
      dplyr::mutate(nvars = dplyr::n()) |>
      dplyr::group_by(target, tf) |>
      dplyr::mutate(tf_sites_per_gene = dplyr::n()) |>
      dplyr::group_by(target) |>
      dplyr::mutate(
        tf_per_gene = length(unique(tf)),
        peak_per_gene = length(unique(region))
      ) |>
      dplyr::group_by(tf) |>
      dplyr::mutate(gene_per_tf = length(unique(target))) |>
      dplyr::group_by(target, tf)

    if (fit_method %in% c("srm", "cv.glmnet", "glmnet", "xgb")) {
      modules <- modules |>
        dplyr::reframe(
          coefficient = sum(coefficient),
          n_regions = peak_per_gene,
          n_genes = gene_per_tf,
          n_tfs = tf_per_gene,
          regions = paste(region, collapse = ";")
        )
    } else {
      modules <- modules |>
        dplyr::reframe(
          coefficient = sum(coefficient),
          n_regions = peak_per_gene,
          n_genes = gene_per_tf,
          n_tfs = tf_per_gene,
          regions = paste(region, collapse = ";"),
          pval = min(pval),
          padj = min(padj)
        )
    }

    modules <- modules |>
      dplyr::distinct() |>
      dplyr::arrange(tf)

    module_pos <- modules |>
      dplyr::filter(coefficient > 0) |>
      dplyr::group_by(tf) |>
      dplyr::filter(dplyr::n() > min_genes_per_module) |>
      dplyr::group_split()
    names(module_pos) <- purrr::map_chr(
      module_pos, function(x) x$tf[[1]]
    )
    module_pos <- purrr::map(module_pos, function(x) x$target)

    module_neg <- modules |>
      dplyr::filter(coefficient < 0) |>
      dplyr::group_by(tf) |>
      dplyr::filter(dplyr::n() > min_genes_per_module) |>
      dplyr::group_split()
    names(module_neg) <- purrr::map_chr(
      module_neg, function(x) x$tf[[1]]
    )
    module_neg <- purrr::map(module_neg, function(x) x$target)

    regions_pos <- modules |>
      dplyr::filter(coefficient > 0) |>
      dplyr::group_by(tf) |>
      dplyr::filter(dplyr::n() > min_genes_per_module) |>
      dplyr::group_split()
    names(regions_pos) <- purrr::map_chr(
      regions_pos, function(x) x$tf[[1]]
    )
    regions_pos <- purrr::map(
      regions_pos, function(x) {
        unlist(stringr::str_split(x$regions, ";"))
      }
    )

    regions_neg <- modules |>
      dplyr::filter(coefficient < 0) |>
      dplyr::group_by(tf) |>
      dplyr::filter(dplyr::n() > min_genes_per_module) |>
      dplyr::group_split()
    names(regions_neg) <- purrr::map_chr(
      regions_neg, function(x) x$tf[[1]]
    )
    regions_neg <- purrr::map(
      regions_neg, function(x) {
        unlist(stringr::str_split(x$regions, ";"))
      }
    )

    module_feats <- list(
      "genes_pos" = module_pos,
      "genes_neg" = module_neg,
      "regions_pos" = regions_pos,
      "regions_neg" = regions_neg
    )

    thisutils::log_message(
      "Found {.val {length(unique(modules$tf))}} TF modules",
      verbose = verbose,
      message_type = "success"
    )

    module_meta <- dplyr::select(
      modules,
      tf,
      target,
      tidyselect::everything()
    )
    object@modules@meta <- module_meta
    object@modules@features <- module_feats
    object@modules@params <- list(
      p_thresh = p_thresh,
      rsq_thresh = rsq_thresh,
      nvar_thresh = nvar_thresh,
      min_genes_per_module = min_genes_per_module
    )

    return(object)
  }
)

#' @param object An object.
#' @param network Name of the network to use.
#'
#' @rdname find_modules
#' @method find_modules CSNObject
#'
#' @export
#' @return A CSNObject object
setMethod(
  f = "find_modules",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        p_thresh = 0.05,
                        rsq_thresh = 0.1,
                        nvar_thresh = 10,
                        min_genes_per_module = 5,
                        ...) {
    params <- Params(object)
    regions <- NetworkRegions(object)
    cell_types <- object@metadata$celltypes

    net_obj <- GetNetwork(
      object,
      network = network,
      celltypes = cell_types
    )

    net_obj <- purrr::map(
      net_obj, function(x) {
        find_modules(
          x,
          p_thresh = p_thresh,
          rsq_thresh = rsq_thresh,
          nvar_thresh = nvar_thresh,
          min_genes_per_module = min_genes_per_module
        )
      }
    )

    modules_list <- purrr::map(
      net_obj, function(x) {
        modules <- x@modules
        reg2peaks <- rownames(
          GetAssay(object, assay = params$peak_assay)
        )[regions@peaks]
        names(reg2peaks) <- Signac::GRangesToString(regions@ranges)
        peaks_pos <- modules@features$regions_pos |>
          purrr::map(
            function(x) unique(reg2peaks[x])
          )
        peaks_neg <- modules@features$regions_neg |>
          purrr::map(
            function(x) unique(reg2peaks[x])
          )
        modules@features[["peaks_pos"]] <- peaks_pos
        modules@features[["peaks_neg"]] <- peaks_neg
        modules
      }
    )

    names(modules_list) <- cell_types

    for (celltype in cell_types) {
      object@networks[[network]][[celltype]]@modules <- modules_list[[celltype]]
    }

    return(object)
  }
)
