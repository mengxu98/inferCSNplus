#' @title Format network coefficients
#'
#' @param coefs A data frame with coefficients
#' @param term term
#' @param adjust_method adjust_method
#'
#' @return A data frame.
#'
#' @export
format_coefs <- function(
    coefs,
    term = ":",
    adjust_method = "fdr") {
  if (dim(coefs)[1] == 0) {
    return(coefs)
  }

  if ("pval" %in% colnames(coefs)) {
    coefs$padj <- stats::p.adjust(
      coefs$pval,
      method = adjust_method
    )
  }

  term_pattern <- paste0("(.+)", term, "(.+)")
  region_pattern <- "[\\d\\w]+_\\d+_\\d+"
  coefs_use <- coefs %>%
    dplyr::filter(!term %in% c("(Intercept)", "Intercept")) %>%
    dplyr::mutate(
      tf_ = stringr::str_replace(term, term_pattern, "\\1"),
      region_ = stringr::str_replace(term, term_pattern, "\\2")
    ) %>%
    dplyr::mutate(
      tf = ifelse(stringr::str_detect(tf_, region_pattern), region_, tf_),
      region = ifelse(!stringr::str_detect(tf_, region_pattern), region_, tf_)
    ) %>%
    dplyr::select(-region_, -tf_) %>%
    dplyr::mutate(
      region = stringr::str_replace_all(region, "_", "-"),
      tf = stringr::str_replace_all(tf, "_", "-"),
      target = stringr::str_replace_all(target, "_", "-")
    ) %>%
    dplyr::select(tf, target, region, term, tidyselect::everything())

  return(coefs_use)
}

#' @param p_thresh Float indicating the significance threshold on the adjusted p-value.
#' @param rsq_thresh Float indicating the \eqn{R^2} threshold on the adjusted p-value.
#' @param nvar_thresh Integer indicating the minimum number of variables in the model.
#' @param min_genes_per_module Integer indicating the minimum number of genes in a module.
#' @param xgb_method Method to get modules from xgb models
#' * \code{'tf'} - Choose top targets for each TF.
#' * \code{'target'} - Choose top TFs for each target gene.
#' @param xgb_top Interger indicating how many top targets/TFs to return.
#' @param verbose Print messages.
#'
#' @return A Network object.
#'
#' @rdname find_modules
#' @export
#' @method find_modules Network
find_modules.Network <- function(
    object,
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

  if (!fit_method %in% c("srm", "glm", "cv.glmnet", "glmnet", "brms", "xgb")) {
    stop(paste0('find_modules() is not yet implemented for "', fit_method, '" models'))
  }

  models_use <- gof(object) %>%
    dplyr::filter(rsq > rsq_thresh & nvariables > nvar_thresh) %>%
    dplyr::pull(target) %>%
    unique()

  modules <- coef(object) %>%
    dplyr::filter(target %in% models_use)

  if (fit_method %in% c("srm", "cv.glmnet", "glmnet")) {
    modules <- modules %>%
      dplyr::filter(estimate != 0)
  } else if (fit_method == "xgb") {
    modules <- modules %>%
      dplyr::group_by_at(xgb_method) %>%
      dplyr::top_n(xgb_top, gain) %>%
      dplyr::mutate(estimate = sign(corr) * gain)
  } else {
    modules <- modules %>%
      dplyr::filter(ifelse(is.na(padj), T, padj < p_thresh))
  }

  modules <- modules %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(nvars = dplyr::n()) %>%
    dplyr::group_by(target, tf) %>%
    dplyr::mutate(tf_sites_per_gene = dplyr::n()) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(
      tf_per_gene = length(unique(tf)),
      peak_per_gene = length(unique(region))
    ) %>%
    dplyr::group_by(tf) %>%
    dplyr::mutate(gene_per_tf = length(unique(target))) %>%
    dplyr::group_by(target, tf)

  if (fit_method %in% c("srm", "cv.glmnet", "glmnet", "xgb")) {
    modules <- modules %>%
      dplyr::reframe(
        estimate = sum(estimate),
        n_regions = peak_per_gene,
        n_genes = gene_per_tf,
        n_tfs = tf_per_gene,
        regions = paste(region, collapse = ";")
      )
  } else {
    modules <- modules %>%
      dplyr::reframe(
        estimate = sum(estimate),
        n_regions = peak_per_gene,
        n_genes = gene_per_tf,
        n_tfs = tf_per_gene,
        regions = paste(region, collapse = ";"),
        pval = min(pval),
        padj = min(padj)
      )
  }

  modules <- modules %>%
    dplyr::distinct() %>%
    dplyr::arrange(tf)

  module_pos <- modules %>%
    dplyr::filter(estimate > 0) %>%
    dplyr::group_by(tf) %>%
    dplyr::filter(dplyr::n() > min_genes_per_module) %>%
    dplyr::group_split() %>%
    {
      names(.) <- purrr::map_chr(., function(x) x$tf[[1]])
      .
    } %>%
    purrr::map(function(x) x$target)

  module_neg <- modules %>%
    dplyr::filter(estimate < 0) %>%
    dplyr::group_by(tf) %>%
    dplyr::filter(dplyr::n() > min_genes_per_module) %>%
    dplyr::group_split() %>%
    {
      names(.) <- purrr::map_chr(., function(x) x$tf[[1]])
      .
    } %>%
    purrr::map(function(x) x$target)

  regions_pos <- modules %>%
    dplyr::filter(estimate > 0) %>%
    dplyr::group_by(tf) %>%
    dplyr::filter(dplyr::n() > min_genes_per_module) %>%
    dplyr::group_split() %>%
    {
      names(.) <- purrr::map_chr(., function(x) x$tf[[1]])
      .
    } %>%
    purrr::map(function(x) unlist(stringr::str_split(x$regions, ";")))

  regions_neg <- modules %>%
    dplyr::filter(estimate < 0) %>%
    dplyr::group_by(tf) %>%
    dplyr::filter(dplyr::n() > min_genes_per_module) %>%
    dplyr::group_split() %>%
    {
      names(.) <- purrr::map_chr(., function(x) x$tf[[1]])
      .
    } %>%
    purrr::map(function(x) unlist(stringr::str_split(x$regions, ";")))

  module_feats <- list(
    "genes_pos" = module_pos,
    "genes_neg" = module_neg,
    "regions_pos" = regions_pos,
    "regions_neg" = regions_neg
  )

  log_message(
    paste0("Found ", length(unique(modules$tf)), " TF modules"),
    verbose = verbose
  )

  module_meta <- dplyr::select(modules, tf, target, tidyselect::everything())
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

#' @param object An object.
#' @param network Name of the network to use.
#'
#' @rdname find_modules
#' @method find_modules CSNObject
#'
#' @export
#' @return A CSNObject object.
find_modules.CSNObject <- function(
    object,
    network = DefaultNetwork(object),
    p_thresh = 0.05,
    rsq_thresh = 0.1,
    nvar_thresh = 10,
    min_genes_per_module = 5,
    ...) {
  params <- Params(object)
  regions <- NetworkRegions(object)
  net_obj <- GetNetwork(object, network = network)
  net_obj <- find_modules(
    net_obj,
    p_thresh = p_thresh,
    rsq_thresh = rsq_thresh,
    nvar_thresh = nvar_thresh,
    min_genes_per_module = min_genes_per_module
  )
  modules <- NetworkModules(net_obj)

  reg2peaks <- rownames(
    GetAssay(object, assay = params$peak_assay)
  )[regions@peaks]
  names(reg2peaks) <- Signac::GRangesToString(regions@ranges)
  peaks_pos <- modules@features$regions_pos %>% purrr::map(function(x) unique(reg2peaks[x]))
  peaks_neg <- modules@features$regions_neg %>% purrr::map(function(x) unique(reg2peaks[x]))
  modules@features[["peaks_pos"]] <- peaks_pos
  modules@features[["peaks_neg"]] <- peaks_neg
  object@csn@networks[[network]]@modules <- modules

  return(object)
}


#' Find peaks or regions near gene body or TSS
#'
#' @param peaks A \code{GRanges} object with peak regions.
#' @param genes A \code{GRanges} object with gene coordinates.
#' @param method Character specifying the method to link peak overlapping motif regions to nearby genes.
#' On of 'Signac' or 'GREAT'.
#' @param upstream Integer defining the distance upstream of the gene/TSS to consider.
#' @param downstream Integer defining the distance downstream of the gene/TSS to consider.
#' @param extend Integer defining the distance from the upstream and downstream of the
#' basal regulatory region. Used only when method is 'GREAT'.
#' @param sep Vector of separators to use for genomic string.
#' First element is used to separate chromosome and coordinates,
#' second separator is used to separate start and end coordinates.
#' @param only_tss Logical value. Measure distance from the TSS (\code{TRUE})
#' or from the entire gene body (\code{FALSE}).
#' @param verbose Logical value. Display messages
#'
#' @return A sparse binary Matrix with gene/peak matches.
#'
#' @export
find_peaks_near_genes <- function(
    peaks,
    genes,
    sep = c("-", "-"),
    method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    verbose = TRUE) {
  method <- match.arg(method)

  if (method == "Signac") {
    if (only_tss) {
      genes <- IRanges::resize(x = genes, width = 1, fix = "start")
    }
    genes_extended <- suppressWarnings(
      expr = Signac::Extend(
        genes,
        upstream = upstream,
        downstream = downstream
      )
    )
    overlaps <- IRanges::findOverlaps(
      query = peaks,
      subject = genes_extended,
      type = "any",
      select = "all"
    )
    hit_matrix <- Matrix::sparseMatrix(
      i = S4Vectors::queryHits(overlaps),
      j = S4Vectors::subjectHits(overlaps),
      x = 1,
      dims = c(
        length(peaks),
        length(genes_extended)
      )
    )
    rownames(hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name
  } else if (method == "GREAT") {
    # Read gene annotation (Ensembl v93, GRCh38)
    utils::data(EnsDb.Hsapiens.v93.annot.UCSC.hg38, envir = environment())
    gene_annot_use <- EnsDb.Hsapiens.v93.annot.UCSC.hg38[
      which(EnsDb.Hsapiens.v93.annot.UCSC.hg38$gene_name %in% genes$gene_name),
    ]
    gene_annot_tss <- dplyr::select(
      tibble::as_tibble(gene_annot_use), seqnames,
      "start" = tss, "end" = tss, strand
    )

    # Create GRanges object storing the TSS information
    tss <- GenomicRanges::GRanges(gene_annot_use)

    # Define basal regulatory region (promoter region)
    # as 5 kb upstream + 1 kb downstream of the TSS
    basal_reg <- suppressWarnings(
      expr = Signac::Extend(
        tss,
        upstream = upstream,
        downstream = downstream
      )
    )

    # Step 1 - get peaks overlap with basal regulatory region
    basal_overlaps <- suppressWarnings(
      IRanges::findOverlaps(
        query = peaks,
        subject = basal_reg,
        type = "any",
        select = "all",
        minoverlap = 2
      )
    )

    peak_all <- Signac::GRangesToString(grange = peaks, sep = sep)
    basal_peak_mapped_idx <- S4Vectors::queryHits(basal_overlaps)
    # basal_mapped_peaks <- unique(peak_all[basal_peak_mapped_idx])
    # n1 <- length(basal_mapped_peaks)

    # Step 2: for the peaks not overlapped with basal regulatory regions,
    # check whether they located within gene body of any genes
    # peak_unmapped_idx <- setdiff(seq(length(peak_all)), basal_peak_mapped_idx)
    peak_unmapped_idx <- setdiff(seq_along(peak_all), basal_peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)

    # Create GRanges object storing annotated gene boundary
    gene_bound <- GenomicRanges::GRanges(gene_annot_use)
    body_overlaps <- IRanges::findOverlaps(
      query = peak_unmapped_region,
      subject = gene_bound,
      type = "any",
      select = "all",
      minoverlap = 2
    )
    body_peak_mapped_idx <- peak_unmapped_idx[S4Vectors::queryHits(body_overlaps)]
    # body_mapped_peaks <- unique(peak_all[body_peak_mapped_idx])
    # n2 <- length(body_mapped_peaks)
    peak_mapped_idx <- c(basal_peak_mapped_idx, body_peak_mapped_idx)

    # Step 3: for the peaks not overlapped with basal regulatory regions of any genes,
    # check whether they overlap with extended regulatory region. i.e. +/- 1MB of basal regulatory region
    peak_unmapped_idx <- setdiff(seq_along(peak_all), peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)
    extend_reg <- suppressWarnings(
      expr = Signac::Extend(
        basal_reg,
        upstream = extend,
        downstream = extend
      )
    )

    # Get overlap between unmapped_peak_region and extended regulatory region
    extended_overlaps <- suppressWarnings(
      IRanges::findOverlaps(
        query = peak_unmapped_region,
        subject = extend_reg,
        type = "any",
        select = "all",
        minoverlap = 2
      )
    )
    extended_peak_mapped_idx <- peak_unmapped_idx[S4Vectors::queryHits(extended_overlaps)]
    # extended_mapped_peaks <- unique(peak_all[extended_peak_mapped_idx])
    # n3 <- length(extended_mapped_peaks)

    hit_matrix <- Matrix::sparseMatrix(
      i = c(
        basal_peak_mapped_idx,
        body_peak_mapped_idx,
        extended_peak_mapped_idx
      ),
      j = c(
        S4Vectors::subjectHits(basal_overlaps),
        S4Vectors::subjectHits(body_overlaps),
        S4Vectors::subjectHits(extended_overlaps)
      ),
      x = 1,
      dims = c(length(peaks), length(basal_reg))
    )
    rownames(hit_matrix) <- peak_all
    colnames(hit_matrix) <- c(basal_reg$gene_name)
  }

  return(hit_matrix)
}

#' Copy of the aggregate.Matrix function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#'
#' @param x x
#' @param groupings groupings
#' @param form form
#' @param fun fun
#' @param ... other parameters
#'
#' @return aggregated matrix
fast_aggregate <- function(
    x,
    groupings = NULL,
    form = NULL,
    fun = "sum",
    ...) {
  if (!is(x, "Matrix")) {
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  }
  if (fun == "count") {
    x <- x != 0
  }
  groupings2 <- groupings
  if (!is(groupings2, "data.frame")) {
    groupings2 <- as.data.frame(groupings2)
  }
  groupings2 <- data.frame(lapply(groupings2, as.factor))
  groupings2 <- data.frame(interaction(groupings2, sep = "_"))
  colnames(groupings2) <- "A"
  if (is.null(form)) {
    form <- stats::as.formula("~0+.")
  }
  form <- stats::as.formula(form)
  mapping <- dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- t(mapping) %*% x
  if (fun == "mean") {
    result@x <- result@x / (fast_aggregate(x, groupings2, fun = "count"))@x
  }
  attr(result, "crosswalk") <- grr::extract(
    groupings, match(rownames(result), groupings2$A)
  )
  return(result)
}

#' Copy of the dMcast function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#'
#' @param data data
#' @param formula formula
#' @param fun.aggregate fun.aggregate
#' @param value.var value.var
#' @param as.factors as.factors
#' @param factor.nas factor.nas
#' @param drop.unused.levels drop.unused.levels
#'
#' @return data
dMcast <- function(
    data,
    formula,
    fun.aggregate = "sum",
    value.var = NULL,
    as.factors = FALSE,
    factor.nas = TRUE,
    drop.unused.levels = TRUE) {
  values <- 1
  if (!is.null(value.var)) {
    values <- data[, value.var]
  }
  alltms <- stats::terms(formula, data = data)
  response <- rownames(attr(alltms, "factors"))[attr(alltms, "response")]
  tm <- attr(alltms, "term.labels")
  interactions <- tm[grep(":", tm)]
  simple <- setdiff(tm, interactions)
  i2 <- strsplit(interactions, ":")
  newterms <- unlist(
    lapply(
      i2, function(x) {
        paste("paste(", paste(x, collapse = ","), ",", "sep='_'", ")")
      }
    )
  )
  newterms <- c(simple, newterms)
  newformula <- stats::as.formula(
    paste("~0+", paste(newterms, collapse = "+"))
  )
  allvars <- all.vars(alltms)
  data <- data[, c(allvars), drop = FALSE]
  if (as.factors) {
    data <- data.frame(lapply(data, as.factor))
  }
  characters <- unlist(lapply(data, is.character))
  data[, characters] <- lapply(data[, characters, drop = FALSE], as.factor)
  factors <- unlist(lapply(data, is.factor))
  # Prevents errors with 1 or fewer distinct levels
  data[, factors] <- lapply(
    data[, factors, drop = FALSE], function(x) {
      if (factor.nas) {
        if (any(is.na(x))) {
          levels(x) <- c(levels(x), "NA")
          x[is.na(x)] <- "NA"
        }
      }
      if (drop.unused.levels) {
        if (nlevels(x) != length(na.omit(unique(x)))) {
          x <- factor(as.character(x))
        }
      }
      y <- stats::contrasts(x, contrasts = FALSE, sparse = TRUE)
      attr(x, "contrasts") <- y
      return(x)
    }
  )
  # Allows NAs to pass
  attr(data, "na.action") <- na.pass
  result <- Matrix::sparse.model.matrix(
    newformula,
    data,
    drop.unused.levels = FALSE,
    row.names = FALSE
  )
  broken_names <- grep("paste(", colnames(result), fixed = TRUE)
  colnames(result)[broken_names] <- lapply(
    colnames(result)[broken_names], function(x) {
      x <- gsub("paste(", replacement = "", x = x, fixed = TRUE)
      x <- gsub(pattern = ", ", replacement = "_", x = x, fixed = TRUE)
      x <- gsub(
        pattern = '_sep = \"_\")',
        replacement = "",
        x = x,
        fixed = TRUE
      )
      return(x)
    }
  )

  result <- result * values
  if (isTRUE(response > 0)) {
    responses <- all.vars(
      stats::terms(stats::as.formula(paste(response, "~0")))
    )
    result <- fast_aggregate(
      result,
      data[, responses, drop = FALSE],
      fun = fun.aggregate
    )
  }
  return(result)
}

#' Aggregate matrix over groups
#'
#' @param x A matrix.
#' @param groups A character vector with the groups to aggregate over.
#' @param fun The summary function to be applied to each group.
#'
#' @return A summary matrix.
#'
#' @export
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = "mean") {
  if (length(groups) == nrow(x) && "character" %in% class(fun)) {
    if (fun %in% c("count", "sum")) {
      agg_mat <- fast_aggregate(x = x, groupings = groups, fun = fun)
      return(agg_mat)
    }

    if (fun == "mean") {
      group_counts <- as.numeric(table(groups))
      agg_mat <- fast_aggregate(x = x, groupings = groups, fun = "sum")
      agg_mat <- agg_mat / group_counts
      return(agg_mat)
    }
  }

  if ("character" %in% class(fun)) {
    fun <- summary_fun[[fun]]
  }

  if (length(groups) == nrow(x)) {
    agg_mat <- sapply(levels(factor(groups)), function(g) {
      chunk <- x[which(groups == g), ]
      if (is.null(dim(chunk))) {
        return(chunk)
      } else {
        return(fun(chunk))
      }
    })
    agg_mat <- Matrix::Matrix(agg_mat, sparse = TRUE)
  } else if (length(groups) <= 1) {
    agg_mat <- fun(x)
    agg_mat <- Matrix::Matrix(agg_mat, sparse = TRUE)
    colnames(agg_mat) <- groups
    rownames(agg_mat) <- colnames(x)
  } else {
    stop("Length of groups must be either nrow(x) or 1.")
  }

  return(Matrix::t(agg_mat))
}


#' Aggregate Seurat assay over groups
#'
#' @param object The csn object.
#' @param group_name A character vector indicating the metadata column to aggregate over.
#' @param fun The summary function to be applied to each group.
#' @param assay The assay to summarize.
#' @param slot The slot to summarize.
#'
#' @return A Seurat object.
#'
#' @export
aggregate_assay <- function(
    object,
    group_name,
    fun = "mean",
    assay = "RNA",
    slot = "data") {
  ass_mat <- Matrix::t(
    Seurat::GetAssayData(object, assay = assay, slot = slot)
  )
  groups <- as.character(object@meta.data[[group_name]])
  agg_mat <- aggregate_matrix(ass_mat, groups = groups, fun = fun)
  if (is.null(object@assays[[assay]]@misc$summary)) {
    object@assays[[assay]]@misc$summary <- list()
  }
  object@assays[[assay]]@misc$summary[[group_name]] <- agg_mat
  return(object)
}
