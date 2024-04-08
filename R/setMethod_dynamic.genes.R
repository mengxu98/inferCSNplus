#' Select dynamic genes
#'
#' @param object Expression object.
#' @param pseudotime Pseudotime
#' @param fdr_threshold Threshold of fdr.
#'
#' @return Gene list
#' @export
#'
#' @method dynamic.genes default
#'
#' @rdname dynamic.genes

#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' vector_result <- inferVECTOR(example_matrix)
#' dynamic.genes(
#'   object = t(vector_result@matrix),
#'   pseudotime = vector_result@pseudotime[, 2]
#' )
dynamic.genes.default <- function(
    object,
    pseudotime = NULL,
    fdr_threshold = 0.05,
    ...) {
  if (is.null(pseudotime)) {
    return(rownames(object))
  }
  sorted_genes <- names(
    sort(apply(object, 1, var), decreasing = TRUE)
  )
  object <- object[sorted_genes, ]

  # Fit a GAM model with a loess term for pseudotime
  res <- purrr::map_dfr(
    sorted_genes,
    function(x) {
      data <- data.frame(
        exp = object[x, ],
        t = pseudotime
      )

      gam_model <- suppressWarnings(
        gam::gam(exp ~ gam::lo(t), data = data)
      )
      data.frame(
        gene = x,
        P_value = summary(gam_model)[4][[1]][1, 5]
      )
    }
  )
  res$fdr <- p.adjust(res$P_value, method = "BH")
  genes <- res[res$fdr < fdr_threshold, 1]
  if (length(genes) < 3) {
    genes <- sorted_genes
  }
  return(genes)
}

gamFit <- function(
    matrix,
    celltime) {
  message("Starting gammma.")
  # could print out if any missing
  ans <- apply(matrix, 1, function(z) {
    d <- data.frame(z = z, t = celltime)
    tmp <- gam::gam(z ~ gam::lo(celltime), data = d)
    p <- summary(tmp)[4][[1]][1, 5]
    p
  })
  ans
}

#' find genes expressed dynamically
#'
#' @param object properly normalized expression matrix
#' @param meta_data sample table that includes pseudotime, rownames = cells, and a group column
#' @param cluster_by vector of group names to include
#' @param group_column column name in meta_data annotating groups in the cluster_by
#' @param pseudotime_column column name in meta_data annotating pseudotime or latent time
#'
#' @return pvals and cell info
#'
#' @export
findDynGenes <- function(
    object,
    meta_data,
    cluster_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime") {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]


  meta_data$cluster <- meta_data[, group_column]
  meta_data$cells <- rownames(meta_data)

  if (is.null(cluster_by)) {
    cluster_by <- unique(meta_data[, group_column])
  }

  meta_data <- purrr::map_dfr(
    cluster_by, function(x) {
      filter(meta_data, cluster == x)
    }
  )

  if (nrow(meta_data) == 0) {
    ans <- list(genes = NULL, cells = NULL)
    return(ans)
  }

  object <- object[, meta_data$cells]

  dynamic_genes <- gamFit(object, meta_data$pseudotime)
  cells <- data.frame(
    cells = meta_data$cells,
    pseudotime = meta_data$pseudotime,
    group = as.vector(meta_data$cluster)
  )
  cells <- cells[order(cells$pseudotime), ]

  ans <- list(genes = dynamic_genes, cells = cells)

  ans
}

#' @export
#'
#' @method dynamic.genes VECTOR
#'
#' @rdname dynamic.genes
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' vector_result <- inferVECTOR(example_matrix)
#' dynamic.genes(vector_result)
dynamic.genes.VECTOR <- function(
    object,
    fdr_threshold = 0.05,
    ...) {
  genes <- dynamic.genes(
    t(object@matrix),
    pseudotime = object@pseudotime[, 2],
    fdr_threshold = fdr_threshold
  )

  return(genes)
}

#' @param cluster cluster
#' @param group_by idents
#' @param slot slot used
#'
#' @export
#'
#' @method dynamic.genes Seurat
#'
#' @rdname dynamic.genes
dynamic.genes.Seurat <- function(
    object,
    fdr_threshold = 0.05,
    group_by = NULL,
    cluster = NULL,
    slot = "counts",
    ...) {
  if (!is.null(group_by)) {
    Seurat::Idents(object) <- group_by
  }
  if (!is.null(cluster)) {
    object$cluster <- Seurat::Idents(object)
    object <- subset(object, cluster == cluster)
  }

  genes <- dynamic.genes(
    Matrix::as.matrix(
      switch(
        EXPR = slot,
        "counts" = object@assays$RNA$counts,
        "data" = object@assays$RNA$data
      )
    ),
    pseudotime = object@meta.data$pseudotime,
    fdr_threshold = fdr_threshold
  )

  return(genes)
}

#' find genes expressed dynamically
#'
#' @param object properly normalized expression matrix
#' @param meta_data sample table that includes pseudotime, rownames = cells, and a group column
#' @param cluster_by vector of group names to include
#' @param group_column column name in meta_data annotating groups in the cluster_by
#' @param pseudotime_column column name in meta_data annotating pseudotime or latent time
#' @param min_cells min cells
#' @param p_value p_value = 0.05
#'
#' @return pvals and cell info
#'
#' @export
findDynGenes_new <- function(
    object,
    meta_data,
    cluster_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime",
    min_cells = 100,
    p_value = 0.05) {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]

  meta_data$cluster <- meta_data[, group_column]
  meta_data$cells <- rownames(meta_data)

  # if (is.null(cluster_by)) { # 这里是用不同细胞类型，后续添加一个参数，不区分细胞类型，使用所有细胞
  #   cluster_by <- unique(meta_data[, group_column])
  # }
  #
  # meta_data_list <- purrr::map(
  #   cluster_by, function(x) {
  #     res <- filter(meta_data, cluster == x)
  #     if (nrow(res) < min_cells) {
  #       return(NULL)
  #     }
  #     return(res)
  #   }
  # )
  # names(meta_data_list) <- cluster_by
  meta_data_list <- dynamic.windowing(meta_data)
  cluster_by <- names(meta_data_list)

  dynamic_genes_list <- lapply(
    meta_data_list, function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      if (nrow(x) < min_cells) {
        return(NULL)
      }

      object <- object[, x$cells]

      dynamic_genes <- gamFit(object, x$pseudotime)
      dynamic_genes <- as.data.frame(dynamic_genes)
      dynamic_genes$genes <- rownames(dynamic_genes)
      dynamic_genes <- dynamic_genes[dynamic_genes$dynamic_genes < p_value, ]
      dynamic_genes <- na.omit(dynamic_genes)

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
  names(dynamic_genes_list) <- cluster_by

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
  if (class(dynamic_object[[1]]) != "list") {
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
        t <- t.test(chunk[gene, ], background[gene, ])
        ans <- data.frame(gene = gene, mean_diff = (t$estimate[1] - t$estimate[2]), pval = t$p.value)
        diffres <- rbind(diffres, ans)
      }

      diffres$padj <- p.adjust(diffres$pval, method = "BH")
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
assign_epochs_new <- function(
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
        t <- t.test(chunk[gene, ], background[gene, ])
        ans <- data.frame(gene = gene, mean_diff = (t$estimate[1] - t$estimate[2]), pval = t$p.value)
        diffres <- rbind(diffres, ans)
      }

      diffres$padj <- p.adjust(diffres$pval, method = "BH")
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
