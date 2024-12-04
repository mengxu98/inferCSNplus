# Functions to define, assign, and extract the dynamic network 

define.epochs <- function(
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
  #   dynamic_object <- list(dynamic_object)
  #   matrix <- list(matrix)
  # }
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

    new_dynRes[[path]] <- path_dyn
  }

  if (length(new_dynRes) == 1) {
    new_dynRes <- new_dynRes[[1]]
  }
  new_dynRes
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
define_epochs <- function(
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
  #   dynamic_object <- list(dynamic_object)
  #   matrix <- list(matrix)
  # }
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

#' Splits data into epochs
#'
#' Splits data into epochs by assigning cells to epochs
#'
#' @param dynamic_object result of running findDynGenes or compileDynGenes
#' @param cuts vector of pseudotime cutoffs
#' @param epoch_names names of resulting epochs, must have length of length(cuts)+1
#'
#' @return updated dynamic_object with epoch column included in dynamic_object$cells
#' @export
split_epochs_by_pseudotime <- function(dynamic_object, cuts, epoch_names = NULL) {
  sampTab <- dynamic_object$cells

  if (max(cuts) > max(sampTab$pseudotime)) {
    stop("Cuts must be within pseudotime.")
  }

  if (!is.null(epoch_names) & (length(epoch_names) != length(cuts) + 1)) {
    stop("Length of epoch_names must be equal to 1+length(cuts).")
  }

  if (is.null(epoch_names)) {
    epoch_names <- paste0("epoch", (1:(length(cuts) + 1)))
  }

  cuts <- c(-0.1, cuts, max(sampTab$pseudotime))
  sampTab$epoch <- NA
  for (i in 2:length(cuts)) {
    sampTab$epoch[(cuts[i - 1] < sampTab$pseudotime) & (sampTab$pseudotime <= cuts[i])] <- epoch_names[i - 1]
  }

  dynamic_object$cells <- sampTab
  dynamic_object
}

#' Splits data into epochs manually
#'
#' Splits data into epochs given group assignment
#'
#' @param dynamic_object result of running findDynGenes or compileDynGenes
#' @param assignment a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch
#'
#' @return updated dynamic_object with epoch column included in dynamic_object$cells
#' @export
split_epochs_by_group <- function(dynamic_object, assignment) {
  sampTab <- dynamic_object$cells
  sampTab$epoch <- NA

  for (e in names(assignment)) {
    sampTab$epoch[sampTab$group %in% assignment[[e]]] <- e
  }

  dynamic_object$cells <- sampTab
  dynamic_object
}

#' @title find_cuts_by_similarity
#'
#' @description
#'  Returns cuts to define epochs via sliding window comparison
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object result of running findDynGenes or define_epochs
#' @param winSize number of cells to each side to compare each cell to
#' @param limit_to vector of genes on which to base epoch cuts, for example, limiting to TFs
#' @param p_value pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into epochs
#' @export
find_cuts_by_similarity <- function(
    matrix,
    dynamic_object,
    winSize = 2,
    limit_to = NULL,
    p_value = 0.05) {
  # limit exp to dynamically expressed genes
  matrix <- matrix[names(dynamic_object$genes[dynamic_object$genes < p_value]), ]

  if (!is.null(limit_to)) {
    limit_to <- intersect(limit_to, rownames(matrix))
    matrix <- matrix[limit_to, ]
  }

  # compute PCC between all cells -- just easier this way
  xcorr <- cor(matrix[, rownames(dynamic_object$cells)])
  pShift <- rep(0, nrow(xcorr) - 1)
  end <- nrow(xcorr)

  for (i in 2:end - 1) {
    past <- seq(i - 1:max(1, i - winSize))
    pastandfuture <- seq(i + 1:min(end, i + winSize))
    closer_to <- which.max(c(mean(xcorr[i, past]), mean(xcorr[i, pastandfuture])))
    pShift[i] <- c(0, 1)[closer_to]
  }

  consecutive_diffs <- diff(pShift)
  cuts_index <- which(consecutive_diffs == 1)

  # cat("Cut indicies: ",cuts_index, "\n")
  cat("Cut points: ", dynamic_object$cells[cuts_index, ]$pseudotime, "\n")

  # It is kinda confusing, but we want the pt of the cell just prior to the determined cutpoint
  # But, because of the way that this is indexed, we don't need to adjust anything
  dynamic_object$cells[cuts_index, ]$pseudotime
}


#' Returns cuts to define epochs
#'
#' Returns cuts to define epochs via clustering
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object result of running findDynGenes or define_epochs
#' @param num_epochs the number of epochs
#' @param limit_to vector of genes on which to base epoch cuts, for example, limiting to TFs
#' @param method what clustering method to use, either 'kmeans' or 'hierarchical'
#' @param p_value pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into epochs
#' @export
find_cuts_by_clustering <- function(
    matrix,
    dynamic_object,
    num_epochs,
    limit_to = NULL,
    method = "kmeans",
    p_value = 0.05) {
  matrix <- matrix[names(dynamic_object$genes[dynamic_object$genes < p_value]), ]
  # cat(nrow(matrix),"\n")

  if (!is.null(limit_to)) {
    limit_to <- intersect(limit_to, rownames(matrix))
    matrix <- matrix[limit_to, ]
  }

  if (method == "kmeans") {
    clustering <- stats::kmeans(t(matrix), num_epochs, iter.max = 100)$cluster

    cuts <- c()
    for (cluster in unique(clustering)) {
      cells <- names(clustering)[clustering == cluster]
      cuts <- c(cuts, max(dynamic_object$cells$pseudotime[rownames(dynamic_object$cells) %in% cells]))
    }
    cuts <- sort(cuts)
    cuts <- cuts[1:length(cuts) - 1]
  } else if (method == "hierarchical") {
    clustering <- stats::hclust(stats::dist(t(matrix)), method = "centroid")
    clusterCut <- stats::cutree(clustering, 3)
    cuts <- c()
    for (cluster in unique(clusterCut)) {
      cells <- names(clusterCut)[clusterCut == cluster]
      cuts <- c(cuts, max(dynamic_object$cells$pseudotime[rownames(dynamic_object$cells) %in% cells]))
    }
    cuts <- sort(cuts)
    cuts <- cuts[1:length(cuts) - 1]
  } else {
    stop("method must be either 'kmeans' or 'hierarchical'.")
  }

  cat("Cut points: ", cuts, "\n")
  cuts
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
#'
assign_epochs <- function(
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
  epoch_names <- unique(dynamic_object$cells$epoch)
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
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
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
        ans <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, ans)
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

#' Divides grnDF into epochs, filters interactions between genes not in same or consecutive epochs
#'
#' @param grnDF result of GRN reconstruction
#' @param epochs result of running assign_epochs
#' @param epoch_network dataframe outlining higher level epoch connectivity (i.e. epoch transition network).
#' If NULL, will assume epochs is ordered linear trajectory
#'
#' @return list of GRNs across epochs and transitions
#' @export
epochGRN <- function(
    grnDF,
    epochs,
    epoch_network = NULL) {
  epochs$mean_expression <- NULL
  all_dyngenes <- unique(unlist(epochs, use.names = FALSE))

  # assign epoch_network if NULL
  if (is.null(epoch_network)) {
    epoch_network <- data.frame(from = character(), to = character())
    for (i in 1:(length(names(epochs)) - 1)) {
      df <- data.frame(from = c(names(epochs)[i]), to = c(names(epochs)[i + 1]))
      epoch_network <- rbind(epoch_network, df)
    }
  }

  # add in self loops to epoch_network (i.e. for interactions ocurring inside same epoch)
  epoch_network <- rbind(
    epoch_network,
    data.frame(from = names(epochs), to = names(epochs))
  )

  # build network for each epoch and transition
  # this will inherently filter out interactions not belonging to same or consecutive epochs (i.e. anything not in epoch_network)
  epoch_network$name <- paste(epoch_network[, 1], epoch_network[, 2], sep = "..")
  GRN <- vector("list", nrow(epoch_network))
  names(GRN) <- epoch_network$name

  print(epoch_network)

  for (t in 1:nrow(epoch_network)) {
    from <- as.character(epoch_network[t, 1])
    to <- as.character(epoch_network[t, 2])

    # TF must be active in source epoch
    temp <- grnDF[grnDF$regulator %in% epochs[[from]], ]

    # For transition network
    if (from != to) {
      # target turns on: target is not active in source epoch but is active in target epoch
      # target turns off: target is active in source epoch but is not active in target epoch
      # remove any other interaction (i.e. interactions that are constant -- target on in both epochs or off in both epochs)
      remove_tgs_in_both <- intersect(epochs[[to]], epochs[[from]])
      remove_tgs_in_neither <- intersect(setdiff(all_dyngenes, epochs[[from]]), setdiff(all_dyngenes, epochs[[to]]))

      temp <- temp[!(temp$TG %in% remove_tgs_in_both), ]
      temp <- temp[!(temp$TG %in% remove_tgs_in_neither), ]
    }
    # Else, For epoch network (non-transition network), keep all interactions as long as TF is active

    GRN[[epoch_network[t, "name"]]] <- temp
  }

  GRN
}




# =================== Old function to cluster and order genes 
# useful for quick and simple epoch assignment, take the place of define_epochs and assign_epochs

#' cluster and order genes
#'
#' @param expSmooth expression matrix
#' @param matrix matrix
#' @param dynamic_object result of running findDynGenes
#' @param pThresh pval threshold
#' @param k number of clusters
#' @param method method
#'
#' @return data.frame of dynamically expressed genes, cluster, peakTime, ordered by peaktime
#' @export
caoGenes <- function(
    expSmooth,
    matrix, # for alt peakTimes
    dynamic_object,
    k = 3,
    pThresh = 0.01,
    method = "kmeans") {
  # Find dyn genes
  genes <- dynamic_object$genes
  genes <- names(genes[which(genes < pThresh)])
  if (length(genes) == 0) {
    cat("no dynamic genes\n")
    ans <- data.frame()
  } else {
    value <- t(scale(t(expSmooth[genes, ])))
    cat("A\n")
    # cluster
    genedist <- utils_myDist(value)
    geneTree <- stats::hclust(genedist, "ave")

    if (method == "kmeans") {
      kmeansRes <- stats::kmeans(genedist, centers = k)
      geneMods2 <- as.character(kmeansRes$cluster)
    } else {
      if (method == "dct") {
        geneMods <- cutreeDynamicTree(geneTree, deepSplit = FALSE, minModuleSize = 50)
        geneMods2 <- as.character(geneMods)
      } else {
        if (method == "pam") {
          geneMods <- pam(genedist, k = k, cluster.only = TRUE)
          geneMods2 <- as.character(geneMods)
        } else {
          geneMods <- stats::cutree(geneTree, k = k)
          geneMods2 <- as.character(geneMods)
        }
      }
    }

    names(geneMods2) <- genes
    geneMods <- geneMods2

    # find max peak
    genesPeakTimes <- apply(expSmooth[genes, ], 1, which.max)
    cat("B\n")
    # the matrix won't have been pruned ordered
    matrix <- matrix[, colnames(expSmooth)]
    ### genesPeakTimesAlt = apply(matrix[genes,], 1, which.max)
    genesPeakTimesRaw <- apply(matrix[genes, ], 1, findTop)
    cat("C\n")

    ### peakTime = names(sort(gemeakTime))

    # order the clusters
    clusterOrder <- vector()
    clusterNames <- unique(geneMods)
    meanPeakTimes <- rep(0, length(clusterNames))
    names(meanPeakTimes) <- clusterNames
    for (clu in clusterNames) {
      cat(clu, "\t")
      cgenes <- names(which(geneMods == clu))
      cat(length(cgenes), "\t")
      meanPeakTimes[clu] <- mean(genesPeakTimes[cgenes])
      cat(meanPeakTimes[clu], "\n")
    }

    meanPeakTimes <- sort(meanPeakTimes)
    cluOrdered <- names(meanPeakTimes)
    # make the data.frame

    ans <- data.frame()
    for (clu in cluOrdered) {
      cat(clu, "\n")
      cgenes <- names(which(geneMods == clu))
      ptX <- genesPeakTimes[cgenes]
      ptX_raw <- genesPeakTimesRaw[cgenes]
      genesOrdered <- names(sort(ptX))

      ### tmpAns<-data.frame(gene = genesOrdered, peakTime = ptX, peakTimeRaw = ptX_raw, epoch = clu)

      tmpAns <- data.frame(gene = genesOrdered, peakTime = ptX[genesOrdered], peakTimeRaw = ptX_raw[genesOrdered], epoch = clu)


      rownames(tmpAns) <- genesOrdered
      ans <- rbind(ans, tmpAns)
    }
  }

  # now, re-label the clusters so they make sense
  epochs <- unique(as.vector(ans$epoch))
  newepochs <- as.vector(ans$epoch)
  new_e <- 1
  for (i in epochs) { # this should go in order they appear
    newepochs[which(ans$epoch == i)] <- new_e
    new_e <- new_e + 1
  }

  ans$epoch <- newepochs
  ans
}


#' Assigns genes to epochs just based on which mean is maximal
#'
#' @param matrix expression matrix
#' @param dynamic_object result of running findDynGenes
#' @param num_epochs num_epochs
#' @param pThresh pThresh
#' @param toScale toScale
#' @param key_word key_word
#'
#' @return data.frame of dynamically expressed genes, cluster, peakTime, ordered by peaktime
#' @export
assign_epochs_simple <- function(
    matrix,
    dynamic_object,
    num_epochs = 3,
    pThresh = 0.01,
    toScale = FALSE,
    key_word = "epoch") {
  # limit exp to dynamically expressed genes
  exp <- matrix[names(dynamic_object$genes[dynamic_object$genes < pThresh]), ]
  # scale the data if needed
  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }

  navg <- ceiling(ncol(exp) * 0.05)

  # compute thresholds for each gene
  # for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5%
  # set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom
  thresholds <- data.frame(gene = rownames(exp), thresh = rep(0, length(rownames(exp))))
  rownames(thresholds) <- thresholds$gene
  for (gene in rownames(exp)) {
    profile <- exp[gene, ][order(exp[gene, ], decreasing = FALSE)]
    bottom <- mean(profile[1:navg])
    top <- mean(profile[(length(profile) - navg):length(profile)])

    thresh <- ((top - bottom) * 0.33) + bottom
    thresholds[gene, "thresh"] <- thresh
  }

  # order cells in exp along pseudotime-- cells ordered in dynamic_object
  t1 <- dynamic_object$cells$pseudotime
  names(t1) <- as.vector(dynamic_object$cells$cell_name)
  sort(t1, decreasing = FALSE)
  exp <- exp[, names(t1)]

  mean_expression <- data.frame(gene = character(), epoch = numeric(), mean_expression = numeric())

  # Divide epochs by pseudotime

  # create empty list to contain active genes
  epoch_names <- paste0(rep(key_word, num_epochs), seq(1:num_epochs))
  epochs <- vector("list", length(epoch_names))
  names(epochs) <- epoch_names

  # determine activity based on average expression in each epoch
  ptmax <- max(dynamic_object$cells$pseudotime)
  ptmin <- min(dynamic_object$cells$pseudotime)
  chunk_size <- (ptmax - ptmin) / num_epochs

  cellsEps <- rep("", length(names(t1)))
  names(cellsEps) <- names(t1)

  for (i in 1:length(epochs)) {
    lower_bound <- ptmin + ((i - 1) * chunk_size)
    upper_bound <- ptmin + (i * chunk_size)
    chunk_cells <- rownames(dynamic_object$cells[dynamic_object$cells$pseudotime >= lower_bound & dynamic_object$cells$pseudotime <= upper_bound, ])
    chunk <- exp[, chunk_cells]

    chunk_df <- data.frame(means = rowMeans(chunk))
    chunk_df <- cbind(chunk_df, thresholds)
    chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

    epochs[[i]] <- rownames(chunk_df[chunk_df$active, ])
    genesPeakTimes <- apply(chunk, 1, which.max)
    gpt <- as.vector(dynamic_object$cells[chunk_cells, ][genesPeakTimes, ]$pseudotime)

    mean_expression <- rbind(
      mean_expression,
      data.frame(
        gene = rownames(chunk),
        epoch = rep(i, length(rownames(chunk))), mean_expression = rowMeans(chunk),
        peakTime = gpt
      )
    )
    cellsEps[chunk_cells] <- epoch_names[i]
  }

  # assign genes to epochs
  genes <- unique(as.vector(mean_expression$gene))
  cat("n genes: ", length(genes), "\n")
  eps <- rep("", length(genes))
  geneEpPT <- rep(0, length(genes))
  epMean <- rep(0, length(genes))

  names(eps) <- genes
  names(geneEpPT) <- genes
  names(epMean) <- genes
  for (gene in genes) {
    x <- mean_expression[mean_expression$gene == gene, ]
    xi <- which.max(x$mean_expression)
    eps[[gene]] <- as.vector(x[xi, ]$epoch)
    geneEpPT[[gene]] <- as.vector(x[xi, ]$peakTime)
    epMean[[gene]] <- max(x$mean_expression)
  }

  geneDF <- data.frame(
    gene = genes,
    epoch = eps,
    peakTime = geneEpPT,
    epMean = epMean,
    pval = dynamic_object$genes[genes]
  )
  cells2 <- dynamic_object$cells[names(t1), ]
  cells2 <- cbind(cells2, epoch = cellsEps)

  list(genes = geneDF, cells = cells2)
}
