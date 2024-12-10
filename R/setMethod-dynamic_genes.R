#' @title Get dynamic genes
#'
#' @param object The input data, a matrix with cells/samples by genes/features or a seurat object.
#' @param ... Arguments for other methods
#'
#' @return A new object with dynamic genes
#' @export
#'
#' @rdname dynamic_genes
setGeneric(
  "dynamic_genes",
  signature = "object",
  function(object, ...) {
    standardGeneric("dynamic_genes")
  }
)

gam_fit <- function(
    matrix,
    pseudotime,
    cores = 1,
    verbose = TRUE,
    adjust_method = "BH") {
  adjust_method <- match.arg(adjust_method, stats::p.adjust.methods)

  res <- parallelize_fun(
    colnames(matrix),
    function(x) {
      gam_model <- suppressWarnings(
        gam::gam(
          exp ~ gam::lo(t),
          data = data.frame(
            exp = matrix[, x],
            t = pseudotime
          )
        )
      )
      data.frame(
        gene = x,
        p_value = summary(gam_model)[4][[1]][1, 5]
      )
    },
    cores = cores,
    verbose = verbose
  ) |> purrr::list_rbind()
  res$adjust_p_value <- stats::p.adjust(
    res$p_value,
    method = adjust_method
  ) |> stats::na.omit()
  res <- res[order(
    as.numeric(res$adjust_p_value),
    decreasing = FALSE
  ), ]

  return(res)
}

#' @param pseudotime Pseudotime of cells, the length and order should be the same as the expression matrix columns.
#' @param adjust_method The method used to calculate adjust P-value.
#' @inheritParams inferCSN
#'
#' @export
#'
#' @rdname dynamic_genes
#'
#' @examples
#' data("example_matrix")
#' data("example_meta_data")
#' dynamic_genes(
#'   object = example_matrix,
#'   pseudotime = example_meta_data$pseudotime
#' )
#'
#' \dontrun{
#' vector_result <- infer_vector(example_matrix)
#' dynamic_genes(
#'   object = vector_result$matrix,
#'   pseudotime = vector_result$pseudotime[, 2]
#' )
#' }
setMethod(
  "dynamic_genes",
  signature = "matrix",
  function(
      object,
      pseudotime = NULL,
      cores = 1,
      verbose = TRUE,
      adjust_method = "BH",
      ...) {
    if (is.null(pseudotime)) {
      log_message(
        "No pseudotime provided, using all genes.",
        verbose = verbose
      )
      return(colnames(object))
    }

    res <- gam_fit(
      object,
      pseudotime,
      cores = cores,
      verbose = verbose,
      adjust_method = adjust_method
    )

    return(res)
  }
)

#' @param pseudotime_column The column name in meta_data annotating pseudotime or latent time.
#' @param layer The layer used in Seurat object.
#'
#' @export
#'
#' @rdname dynamic_genes
setMethod(
  "dynamic_genes",
  signature = "Seurat",
  function(object,
           pseudotime_column,
           cores = 1,
           verbose = TRUE,
           adjust_method = "BH",
           layer = "data",
           ...) {
    if (!pseudotime_column %in% colnames(object@meta.data)) {
      stop("pseudotime_column not existed in the provides Seurat object.")
    }

    Seurat::Misc(object, "dynamic_genes") <- dynamic_genes(
      Matrix::as.matrix(
        Seurat::GetAssayData(object, layer = layer)
      ),
      pseudotime = object@meta.data[[pseudotime_column]],
      cores = cores,
      verbose = verbose,
      adjust_method = adjust_method,
      ...
    )

    return(object)
  }
)

#' find genes expressed dynamically
#'
#' @param object properly normalized expression matrix
#' @param meta_data sample table that includes pseudotime, rownames = cells, and a group column
#' @param celltype_by vector of group names to include
#' @param group_column column name in meta_data annotating groups in the celltype_by
#' @param pseudotime_column column name in meta_data annotating pseudotime or latent time
#' @param cores cores
#'
#' @return pvals and cell info
#'
#' @export
findDynGenes <- function(
    object,
    meta_data,
    celltype_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime",
    cores = 1) {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]

  meta_data$cluster <- meta_data[, group_column]
  meta_data$cells <- rownames(meta_data)

  if (is.null(celltype_by)) {
    celltype_by <- unique(meta_data[, group_column])
  }

  meta_data <- purrr::map_dfr(
    celltype_by, function(x) {
      filter(meta_data, cluster == x)
    }
  )

  if (nrow(meta_data) == 0) {
    res <- list(genes = NULL, cells = NULL)
    return(res)
  }

  object <- object[, meta_data$cells]

  res <- gam_fit(
    object,
    meta_data$pseudotime,
    cores = cores
  )
  p_value <- res$adjust_p_value
  names(p_value) <- res$gene

  cells <- data.frame(
    cells = meta_data$cells,
    pseudotime = meta_data$pseudotime,
    group = as.vector(meta_data$cluster)
  )
  cells <- cells[order(cells$pseudotime), ]

  return(
    list(
      genes = p_value,
      cells = cells
    )
  )
}

#' find genes expressed dynamically
#'
#' @param object properly normalized expression matrix
#' @param meta_data sample table that includes pseudotime, rownames = cells, and a group column
#' @param celltype_by vector of group names to include
#' @param group_column column name in meta_data annotating groups in the celltype_by
#' @param pseudotime_column column name in meta_data annotating pseudotime or latent time
#' @param min_cells min cells
#' @param p_value p_value = 0.05
#' @param verbose verbose
#' @param cores cores
#'
#' @return pvals and cell info
#'
#' @export
dynamic_genes_new <- function(
    object,
    meta_data,
    celltype_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime",
    min_cells = 100,
    p_value = 0.05,
    verbose = TRUE,
    cores = 1) {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]

  meta_data$cluster <- meta_data[, group_column]
  meta_data$cells <- rownames(meta_data)

  meta_data_list <- dynamic_windowing(
    meta_data,
    group_column = group_column,
    pseudotime_column = pseudotime_column,
    min_cells = min_cells
  )

  celltype_by <- names(meta_data_list)
  dynamic_genes_list <- purrr::map(
    celltype_by, function(x) {
      cluster <- x
      x <- meta_data_list[[x]]
      if (is.null(x)) {
        log_message(cluster, " not existed.", verbose = verbose)
        return(NULL)
      }

      object <- object[, x$cells]

      log_message("\rStarting gammma for cluster: ", cluster, verbose = verbose)
      res <- gam_fit(
        object,
        x$pseudotime,
        cores = cores
      )
      dynamic_genes <- res$adjust_p_value
      names(dynamic_genes) <- res$gene
      dynamic_genes$genes <- rownames(dynamic_genes)
      dynamic_genes <- dynamic_genes[dynamic_genes$dynamic_genes < p_value, ]

      cells <- data.frame(
        cells = x$cells,
        pseudotime = x$pseudotime,
        group = as.vector(x$cluster)
      )
      cells <- cells[order(cells$pseudotime), ]
      res <- list(genes = dynamic_genes, cells = cells)

      return(res)
    }
  )

  names(dynamic_genes_list) <- celltype_by

  meta_data_list <- meta_data_list[!purrr::map_lgl(meta_data_list, is.null)]
  dynamic_genes_list <- dynamic_genes_list[!purrr::map_lgl(dynamic_genes_list, is.null)]

  result <- purrr::map2(
    meta_data_list, dynamic_genes_list, function(x, y) {
      list(x, y)
    }
  )

  return(result)
}

#' Define epochs
#'
#' @param dynamic_object results of running findDynGenes, or a list of results of running findDynGenes per path. If list, names should match names of matrix.
#' @param matrix genes-by-cells expression matrix, or a list of expression matrices per path. If list, names should match names of dynamic_object.
#' @param method method to define epochs. Either "pseudotime", "cell_order", "group", "con_similarity", "kmeans", "hierarchical"
#' @param num_epochs number of epochs to define. Ignored if epoch_transitions, pseudotime_cuts, or group_assignments are provided.
#' @param pseudotime_cuts vector of pseudotime cutoffs. If NULL, cuts are set to max(pseudotime)/num_epochs.
#' @param group_assignments a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch
#' @param p_value p_value
#' @param winSize winSize
#'
#' @return updated list of dynamic_object with epoch column included in dynamic_object$cells
#' @export
define_epochs_new <- function(
    dynamic_object,
    matrix,
    method = "pseudotime",
    num_epochs = 2,
    pseudotime_cuts = NULL,
    group_assignments = NULL,
    p_value = 0.05,
    winSize = 2) {
  # put dynamic_object and matrix into list, if not already
  # for cleaner code later on... just put everything is in list.
  # if (class(dynamic_object[[1]]) != "list") {
  if (!is.list(dynamic_object[[1]])) {
    dynamic_object <- list(dynamic_object)
    matrix <- list(matrix)
  }

  new_dynRes <- list()
  for (path in 1:length(dynamic_object)) {
    if (!is.null(names(dynamic_object))) {
      path <- names(dynamic_object)[path]
    }
    path_dyn <- dynamic_object[[path]]
    path_dyn$cells$epoch <- NA

    # Define epochs by pseudotime: use pseudotime cutoffs
    if (method == "pseudotime") {
      if (is.null(pseudotime_cuts)) {
        # if pseudotime_cuts NULL, then split pseudotime evenly (note: this is not the same as using cell_order)
        pseudotime_cuts <- seq(
          min(path_dyn$cells$pseudotime),
          max(path_dyn$cells$pseudotime),
          (max(path_dyn$cells$pseudotime) - min(path_dyn$cells$pseudotime)) / num_epochs
        )
        pseudotime_cuts <- pseudotime_cuts[-length(pseudotime_cuts)]
        pseudotime_cuts <- pseudotime_cuts[-1]
      }

      path_dyn <- split_epochs_by_pseudotime(path_dyn, pseudotime_cuts)
    }

    # Define epochs by cell order: equal number of cells in each epoch
    if (method == "cell_order") {
      t1 <- path_dyn$cells$pseudotime
      names(t1) <- as.vector(path_dyn$cells$cell_name)
      # t1 <- sort(t1, decreasing = FALSE)
      chunk_size <- floor(length(t1) / num_epochs)

      for (i in 1:num_epochs) {
        if (i == num_epochs) {
          cells_in_epoch <- names(t1)[(1 + ((i - 1) * chunk_size)):length(t1)]
        } else {
          cells_in_epoch <- names(t1)[(1 + ((i - 1) * chunk_size)):(i * chunk_size)]
        }
        path_dyn$cells$epoch[path_dyn$cells$cell_name %in% cells_in_epoch] <- paste0("epoch", i)
      }
    }

    # Define epochs by group: use clusters to separate cells into epochs
    if (method == "group") {
      if (is.null(group_assignments)) {
        stop("Must provide group_assignments for group method.")
      }

      path_dyn <- split_epochs_by_group(path_dyn, group_assignments)
    }

    # define by consecutive similarity
    if (method == "con_similarity") {
      cuts <- find_cuts_by_similarity(matrix[[path]], path_dyn, winSize = winSize, p_value = p_value)
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }

    # define by kmeans
    if (method == "kmeans") {
      cuts <- find_cuts_by_clustering(
        matrix[[path]],
        path_dyn,
        num_epochs = num_epochs,
        method = "kmeans",
        p_value = p_value
      )
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }

    # define by hierarchical clustering
    if (method == "hierarchical") {
      cuts <- find_cuts_by_clustering(matrix[[path]], path_dyn, num_epochs = num_epochs, method = "hierarchical", p_value = p_value)
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }


    new_dynRes[[path]] <- path_dyn
  }

  if (length(new_dynRes) == 1) {
    new_dynRes <- new_dynRes[[1]]
  }
  new_dynRes
}

#' Assigns genes to epochs
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_epochs
#' @param method method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)
#' @param p_value p value threshold if gene is dynamically expressed
#' @param pThresh_DE p value if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.
#'
#' @return epochs a list detailing genes active in each epoch
#' @export
assign_epochs_new1 <- function(
    matrix,
    dynamic_object,
    method = "active_expression",
    p_value = 0.05,
    pThresh_DE = 0.05,
    active_thresh = 0.33,
    toScale = FALSE,
    forceGenes = TRUE) {
  if (active_thresh < 0 | active_thresh > 1) {
    stop("active_thresh must be between 0 and 1.")
  }
  dynamic_object <- dynamic_object[[2]]
  # limit exp to dynamically expressed genes
  exp <- Matrix::as.matrix(
    matrix[intersect(rownames(matrix), names(dynamic_object$genes[dynamic_object$genes < p_value])), ]
  )
  # scale the data if needed
  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }

  # Assign epochs based on assignment in dynamic_object$cells$epoch
  # create empty list to contain active genes
  # epoch_names <- unique(dynamic_object$cells$epoch)
  epoch_names <- unique(dynamic_object$cells$group)
  epochs <- vector("list", length(epoch_names))
  names(epochs) <- epoch_names

  navg <- ceiling(ncol(exp) * 0.05)
  # compute thresholds for each gene
  # for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5%
  # set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom
  thresholds <- data.frame(
    gene = rownames(exp),
    thresh = rep(0, length(rownames(exp)))
  )
  rownames(thresholds) <- thresholds$gene
  for (gene in rownames(exp)) {
    profile <- exp[gene, ][order(exp[gene, ], decreasing = FALSE)]
    bottom <- mean(profile[1:navg])
    top <- mean(profile[(length(profile) - navg):length(profile)])

    thresh <- ((top - bottom) * active_thresh) + bottom
    thresholds[gene, "thresh"] <- thresh
  }

  mean_expression <- data.frame(
    gene = character(),
    epoch = numeric(),
    mean_expression = numeric()
  )
  if (method == "active_expression") {
    # determine activity based on average expression in each epoch
    for (epoch in names(epochs)) {
      # chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == epoch, "cells"]
      chunk <- exp[, chunk_cells]

      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- rownames(chunk_df[chunk_df$active, ])

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  } else {
    # determine which genes are differentially expressed per epoch

    # DESeq will require raw count data
    # for (epoch in names(epochs)){
    # 	st<-dynamic_object$cells
    # 	st$epoch[st$epoch!=epoch]<-"background"
    # 	exp<-exp[,rownames(st)]

    # 	cds<-newCountDataSet(exp,st$epoch)
    # 	cds<-estimateSizeFactors(cds)
    # 	cds<-tryCatch({estimateDispersions(cds)},error=function(e){print("using fitType=local"); estimateDispersions(cds,fitType='local')})
    # 	diffres<-nbinomTest(cds,epoch,"background")

    #   		epochs[[epoch]]<-diffres$id[diffres$padj<pThresh_DE]

    # 	chunk_cells<-rownames(dynamic_object$cells[dynamic_object$cells$epoch==epoch,])
    # 	chunk<-exp[,chunk_cells]
    # 	mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(epoch,length(rownames(chunk))),mean_expression=rowMeans(chunk)))
    # }

    # Data is scaled and log-normalized. Use t-test here
    for (epoch in names(epochs)) {
      # chunk_cells <- rownames(dynamic_object$cells[dynamic_object$cells$epoch == epoch, ])
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
      chunk <- exp[, chunk_cells]

      background <- exp[, !(colnames(exp) %in% chunk_cells)]

      diffres <- data.frame(gene = character(), mean_diff = double(), pval = double())
      for (gene in rownames(exp)) {
        t <- stats::t.test(chunk[gene, ], background[gene, ])
        res <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, res)
      }

      diffres$padj <- stats::p.adjust(diffres$pval, method = "BH")
      diffres <- diffres[diffres$mean_diff > 0, ] # Filter for genes that are on

      epochs[[epoch]] <- diffres$gene[diffres$padj < pThresh_DE]

      # Filter for genes that are actively expressed
      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- intersect(epochs[[epoch]], rownames(chunk_df[chunk_df$active, ]))

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  }

  # some dynamically assigned genes may not be assigned to any epoch
  # if forceGenes, then assign these orphan genes to the epoch in which they have max expression
  if (forceGenes) {
    assignedGenes <- unique(unlist(epochs))
    orphanGenes <- setdiff(rownames(exp), assignedGenes)
    message("There are ", length(orphanGenes), " orphan genes\n")
    for (oGene in orphanGenes) {
      xdat <- mean_expression[mean_expression$gene == oGene, ]
      ep <- xdat[which.max(xdat$mean_expression), ]$epoch
      epochs[[ep]] <- append(epochs[[ep]], oGene)
    }
  }

  epochs$mean_expression <- mean_expression

  epochs
}

#' Assigns genes to epochs
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_epochs
#' @param method method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)
#' @param p_value pval threshold if gene is dynamically expressed
#' @param pThresh_DE pval if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.
#'
#' @return epochs a list detailing genes active in each epoch
#' @export
assign_network <- function(
    matrix,
    dynamic_object,
    method = "active_expression",
    p_value = 0.05,
    pThresh_DE = 0.05,
    active_thresh = 0.33,
    toScale = FALSE,
    forceGenes = TRUE) {
  if (active_thresh < 0 | active_thresh > 1) {
    stop("active_thresh must be between 0 and 1.")
  }
  message("Starting.")

  dynamic_object <- dynamic_object[[2]]
  # limit exp to dynamically expressed genes
  dynamic_genes <- dynamic_object$genes
  dynamic_genes <- dynamic_genes[dynamic_genes$dynamic_genes < p_value, ]$genes
  exp <- Matrix::as.matrix(
    matrix[intersect(rownames(matrix), dynamic_genes), ]
  )
  # scale the data if needed
  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }

  # Assign epochs based on assignment in dynamic_object$cells$epoch
  # create empty list to contain active genes
  epoch_names <- unique(dynamic_object$cells$group)
  epochs <- vector("list", length(epoch_names))
  names(epochs) <- epoch_names

  navg <- ceiling(ncol(exp) * 0.05)
  # compute thresholds for each gene
  # for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5%
  # set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom
  thresholds <- data.frame(
    gene = rownames(exp),
    thresh = rep(0, length(rownames(exp)))
  )
  rownames(thresholds) <- thresholds$gene
  for (gene in rownames(exp)) {
    profile <- exp[gene, ][order(exp[gene, ], decreasing = FALSE)]
    bottom <- mean(profile[1:navg])
    top <- mean(profile[(length(profile) - navg):length(profile)])

    thresh <- ((top - bottom) * active_thresh) + bottom
    thresholds[gene, "thresh"] <- thresh
  }

  mean_expression <- data.frame(
    gene = character(),
    epoch = numeric(),
    mean_expression = numeric()
  )
  if (method == "active_expression") {
    # determine activity based on average expression in each epoch
    for (epoch in names(epochs)) {
      # chunk_cells <- rownames(dynamic_object$cells[dynamic_object$cells$epoch == epoch, ])
      # chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == epoch, "cells"]
      # chunk_cells <- dynamic_object$cells[, "cells"]
      chunk <- exp[, chunk_cells]

      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- rownames(chunk_df[chunk_df$active, ])

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  } else {
    # Data is scaled and log-normalized. Use t-test here
    for (epoch in names(epochs)) {
      # chunk_cells <- rownames(dynamic_object$cells[dynamic_object$cells$epoch == epoch, ])
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
      chunk <- exp[, chunk_cells]

      background <- exp[, !(colnames(exp) %in% chunk_cells)]

      diffres <- data.frame(gene = character(), mean_diff = double(), pval = double())
      for (gene in rownames(exp)) {
        t <- stats::t.test(chunk[gene, ], background[gene, ])
        res <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, res)
      }

      diffres$padj <- stats::p.adjust(diffres$pval, method = "BH")
      diffres <- diffres[diffres$mean_diff > 0, ] # Filter for genes that are on

      epochs[[epoch]] <- diffres$gene[diffres$padj < pThresh_DE]

      # Filter for genes that are actively expressed
      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- intersect(epochs[[epoch]], rownames(chunk_df[chunk_df$active, ]))

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  }

  # some dynamically assigned genes may not be assigned to any epoch
  # if forceGenes, then assign these orphan genes to the epoch in which they have max expression
  if (forceGenes) {
    assignedGenes <- unique(unlist(epochs))
    orphanGenes <- setdiff(rownames(exp), assignedGenes)
    message("There are ", length(orphanGenes), " orphan genes\n")
    for (oGene in orphanGenes) {
      xdat <- mean_expression[mean_expression$gene == oGene, ]
      ep <- xdat[which.max(xdat$mean_expression), ]$epoch
      epochs[[ep]] <- append(epochs[[ep]], oGene)
    }
  }

  epochs$mean_expression <- mean_expression

  epochs
}
