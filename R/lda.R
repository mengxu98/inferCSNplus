construct_lda_data <- function(
    network_list,
    celltypes_order = NULL,
    combine_method = "intersect") {
  networks_name <- names(network_list)
  if (is.null(celltypes_order)) {
    celltypes_order <- networks_name
  } else {
    celltypes_order <- intersect(celltypes_order, networks_name)
  }
  network_list <- network_list[celltypes_order]
  tfs_list <- purrr::map(
    celltypes_order,
    .f = function(x) {
      regulators <- network_list[[x]][, 1]
      regulators[!duplicated(regulators)]
    }
  )

  combine_method <- match.arg(combine_method, c("intersect", "c"))
  if (combine_method == "intersect") {
    tfs_list_intersect <- purrr::reduce(tfs_list, intersect)
    tfs_list <- purrr::map(
      celltypes_order,
      .f = function(x) {
        tfs_list_intersect
      }
    )
  }

  start_index <- c()
  end_index <- c()
  for (i in 1:length(tfs_list)) {
    if (i == 1) {
      start_index <- 1
      end_index <- length(tfs_list[[i]])
    } else {
      start_index <- c(start_index, end_index[i - 1] + 1)
      end_index <- c(end_index, length(purrr::reduce(tfs_list[1:i], c)))
    }
  }

  networks_matrix <- purrr::map2_dfr(
    celltypes_order,
    tfs_list,
    .f = function(x, y) {
      network <- network_list[[x]]
      network_matrix <- table.to.matrix(network, regulators = y)
      as.data.frame(network_matrix[y, ])
    }
  )
  networks_matrix[is.na(networks_matrix)] <- 0
  rownames(networks_matrix) <- NULL

  lad_data <- list(
    networks_matrix = networks_matrix,
    celltypes_order = celltypes_order,
    tfs_list = purrr::reduce(tfs_list, c),
    start_index = start_index,
    end_index = end_index
  )

  return(lad_data)
}

kl_divergence <- function(P, Q) {
  # dist = kl_divergence(P,Q) Kullback-Leibler divergence of two discrete probability
  # distributions
  # P and Q  are automatically normalised to have the sum of one on rows
  # have the length of one at each
  # P =  n x nbins
  # Q =  1 x nbins or n x nbins(one to one)
  # dist = n x 1
  if (ncol(P) != ncol(Q)) {
    stop("The number of columns in P and Q should be the same")
  }

  if (any(!is.finite(P)) || any(!is.finite(Q))) {
    stop("The inputs contain non-finite values!")
  }

  if (nrow(Q) == 1) {
    Q <- Q / sum(Q)
    P <- P / rowSums(P)
    temp <- P * log(P / matrix(Q, nrow = nrow(P), ncol = ncol(P), byrow = TRUE))
    temp[is.na(temp)] <- 0
    dist <- rowSums(temp)
  } else if (nrow(Q) == nrow(P)) {
    Q <- Q / rowSums(Q)
    P <- P / rowSums(P)
    temp <- P * log(P / Q)
    temp[is.na(temp)] <- 0
    dist <- rowSums(temp)
  } else {
    stop("The dimensions of P and Q are not compatible")
  }

  return(dist)
}

lda_analysis <- function(
    lad_data,
    k,
    path = NULL,
    binary = FALSE,
    binary_threshold = 0.5,
    ntop = 30,
    custom_colors = RColorBrewer::brewer.pal(9, "Blues"),
    legend_width = 2,
    legend_font_size = 0.5,
    font_size = 0.7,
    bar_color = "#3366cc",
    line_color = "#3366cc",
    seed = 2024) {
  celltypes <- lda_data$celltypes_order
  networks_matrix <- lda_data$networks_matrix
  start_index <- lda_data$start_index
  end_index <- lda_data$end_index
  regulators <- lda_data$tfs_list
  # idrm <- c()

  # start_index <- c()
  # end_index <- c()
  # regulators <- c()
  # for (i in seq_along(celltypes)) {
  #   cell <- celltypes[i]
  #   dt <- read.delim(
  #     file.path(indir, paste0(cell, "_consensus_edges", prefix, "_mat", cf, ".txt")),
  #     header = TRUE
  #   )
  #
  #   start_index <- c(start_index, (i - 1) * nrow(dt) + 1)
  #   networks_matrix <- rbind(networks_matrix, dt[, -1])
  #   end_index <- c(end_index, nrow(networks_matrix))
  #   regulators <- c(regulators, dt[, 1])
  # }
  networks_matrix <- abs(networks_matrix)
  rowsum <- rowSums(networks_matrix)
  idrms <- which(rowsum == 0)

  # maybe not reasonalbe?
  if (length(idrms) != 0) {
    networks_matrix <- networks_matrix + 0.00001
  }
  rowsum <- rowSums(networks_matrix)
  idrms <- which(rowsum == 0)
  # idrms_f <- which(rowsum != 0)

  # n <- nrow(dt)
  # v <- ncol(dt)
  # regnames <- colnames(dt)[2:length(colnames(dt))]
  # gnames <- dt[, 1]

  set.seed(seed)
  if (is.null(path)) {
    path <- "lda/"
  }
  outdir <- file.path(path, paste0("k", k))
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (binary) {
    # nsd <- 1
    # binary_threshold <- mean(abs(as.matrix(networks_matrix))) + (nsd * sd(abs(as.matrix(networks_matrix))))
    networks_matrix[networks_matrix >= binary_threshold] <- 1
    networks_matrix[networks_matrix < binary_threshold] <- 0
    lda_result <- topicmodels::LDA(networks_matrix, k)
  } else {
    # lda_result <- topicmodels::LDA(ceiling(networks_matrix * 100), k)
    lda_result <- topicmodels::LDA(ceiling(networks_matrix), k)
  }
  # if (binary == 1) {
  #     lda_result <- LDA(networks_matrix[idrms_f,], k)
  # } else {
  #     lda_result <- LDA(ceiling(networks_matrix[idrms_f,] * 100), k)
  # }

  save(
    lda_result,
    file = file.path(outdir, paste0("lda_model_k", k, ".RData"))
  )

  # Plotting Document-Topic and Topic-Word Probabilities
  doc_prob <- modeltools::posterior(lda_result)$topics
  rownames(doc_prob) <- regulators
  word_prob <- modeltools::posterior(lda_result)$terms

  pdf(
    file = file.path(outdir, paste0("lda_k", k, "_allcells", ".pdf")),
    width = 8,
    height = 5
  )
  par(mfrow = c(1, 2))
  image(
    t(doc_prob),
    main = "Document-Topic Probabilities",
    xlab = "Topics",
    ylab = "TFs * cells",
    col = custom_colors,
    axes = FALSE
  )
  axis(
    1,
    at = seq(0, 1, length.out = ncol(doc_prob)),
    labels = colnames(doc_prob),
    cex.axis = font_size
  )
  fields::image.plot(
    add = TRUE,
    legend.only = TRUE,
    zlim = range(doc_prob),
    col = custom_colors,
    legend.width = legend_width,
    legend.cex = legend_font_size
  )
  image(
    word_prob,
    main = "Topic-gene Probabilities",
    xlab = "Topics",
    ylab = "Genes",
    col = custom_colors,
    axes = FALSE
  )
  axis(
    1,
    at = seq(0, 1, length.out = nrow(word_prob)),
    labels = rownames(word_prob),
    cex.axis = font_size
  )
  fields::image.plot(
    add = TRUE,
    legend.only = TRUE,
    zlim = range(doc_prob),
    col = custom_colors,
    legend.width = legend_width,
    legend.cex = legend_font_size
  )
  dev.off()

  docprob <- doc_prob
  topicid <- apply(docprob, 1, which.max)
  tforder <- order(topicid)

  pdf(
    file = file.path(outdir, paste0("lda_k", k, "_tf-topicPercell.pdf")),
    width = 15,
    height = 10
  )
  par(mfrow = c(2, 4))
  for (i in seq_along(celltypes)) {
    id <- start_index[i]:end_index[i]
    idrm <- intersect(idrms, id)
    idkp <- setdiff(id, idrm)
    idkp <- intersect(idkp, seq_len(nrow(docprob)))
    # idkp <- intersect(idrms_f, id)

    d1prob <- docprob[idkp, ]
    tfs <- regulators[idkp]
    tfid <- apply(d1prob, 1, which.max)
    tforder <- order(tfid)
    image(
      t(d1prob[tforder, ]),
      main = celltypes[i],
      xlab = "Topics",
      ylab = "TFs * cells",
      col = custom_colors,
      axes = FALSE
    )
    axis(
      1,
      at = seq(0, 1, length.out = ncol(doc_prob)),
      labels = colnames(doc_prob),
      cex.axis = font_size
    )
    fields::image.plot(
      add = TRUE,
      legend.only = TRUE,
      zlim = range(doc_prob),
      col = custom_colors,
      legend.width = legend_width,
      legend.cex = legend_font_size
    )
  }
  dev.off()

  # Save genes and topic ids
  TopicWordProb <- apply(t(word_prob), 1, max)
  geneid <- apply(t(word_prob), 1, which.max)
  write.table(
    data.frame(names(TopicWordProb), geneid, TopicWordProb),
    file = file.path(outdir, "genes_topicid.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  # Get TF-topic per cell type and save results
  for (i in seq_along(celltypes)) {
    id <- start_index[i]:end_index[i]
    idrm <- intersect(idrms, id)
    idkp <- setdiff(id, idrm)

    TopicProb <- apply(docprob[idkp, ], 1, max)
    TFtopicid <- apply(docprob[idkp, ], 1, which.max)

    tf_topic_file <- file.path(
      outdir,
      paste0("TFs_topicid_", celltypes[i], ".txt")
    )
    write.table(
      data.frame(names(TopicProb), TFtopicid, TopicProb),
      file = tf_topic_file,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE
    )

    d1prob <- docprob[idkp, ]
    tfs <- regulators[idkp]
    tfid <- apply(d1prob, 1, which.max)
    tforder <- order(tfid)

    pdf(
      file = file.path(outdir, paste0("lda_k", k, "_tf-topicPercell.pdf")),
      width = 15,
      height = 10
    )
    par(mfrow = c(2, 4))
    for (i in seq_along(celltypes)) {
      id <- start_index[i]:end_index[i]
      idrm <- intersect(idrms, id)
      idkp <- setdiff(id, idrm)

      d1prob <- docprob[idkp, ]
      tfs <- regulators[idkp]
      tfid <- apply(d1prob, 1, which.max)
      tforder <- order(tfid)

      image(
        t(d1prob[tforder, ]),
        main = celltypes[i],
        xlab = "Topics",
        ylab = "TFs * cells",
        col = custom_colors,
        axes = FALSE
      )
      axis(
        1,
        at = seq(0, 1, length.out = ncol(d1prob[tforder, ])),
        labels = colnames(d1prob[tforder, ]),
        cex.axis = font_size
      )

      fields::image.plot(
        add = TRUE,
        legend.only = TRUE,
        zlim = range(d1prob[tforder, ]),
        col = custom_colors,
        legend.width = legend_width,
        legend.cex = legend_font_size
      )
    }
    dev.off()
  }

  topic_tfs_list <- list()
  genes_select_tables <- as.data.frame(matrix(data = 1, nrow = ntop))
  # Comparing TF-topic distribution between cell types
  for (i in seq_along(celltypes)) {
    id1 <- start_index[i]:end_index[i]
    d1prob <- docprob[id1, ]
    data1_reg <- networks_matrix[id1, ]
    rownames(data1_reg) <- regulators[id1]
    data1_reg_degree <- rowSums(data1_reg)

    for (j in (i + 1):length(celltypes)) {
      if (i == length(celltypes)) {
        j <- 1
      }
      id2 <- start_index[j]:end_index[j]
      d2prob <- docprob[id2, ]
      data2_reg <- networks_matrix[id2, ]
      rownames(data2_reg) <- regulators[id2]
      data2_reg_degree <- rowSums(data2_reg)

      ddist1 <- kl_divergence(d1prob, d2prob)
      ddist2 <- kl_divergence(d2prob, d1prob)
      mydist <- (ddist1 + ddist2) / 2

      rorder <- order(mydist, decreasing = TRUE)
      mydist <- mydist[rorder]

      genes_select <- names(mydist)[1:ntop]

      genes_select_table <- data.frame(
         genes_select
      )
      names(genes_select_table) <- paste0(celltypes[i], "_", celltypes[j])
      genes_select_tables <- cbind(genes_select_tables, genes_select_table)

      pdf(
        file = file.path(
          outdir,
          paste0(
            "lda_k",
            k,
            "_",
            celltypes[i],
            "_",
            celltypes[j],
            "_top",
            ntop,
            "regulator_heatmap_numberoftargets.pdf"
          )
        ),
        width = 18,
        height = 8
      )

      par(mfrow = c(1, 5))

      image(
        t(d1prob[genes_select, ]),
        main = paste(celltypes[i], "regulators"),
        col = custom_colors,
        axes = FALSE
      )
      axis(
        1,
        at = seq(0, 1, length.out = ncol(d1prob)),
        labels = colnames(d1prob),
        cex.axis = font_size
      )
      axis(
        2,
        at = seq(0, 1, length.out = ntop),
        labels = genes_select,
        las = 2,
        cex.axis = font_size
      )

      fields::image.plot(
        add = TRUE,
        legend.only = TRUE,
        zlim = range(d1prob),
        col = custom_colors,
        legend.width = legend_width,
        legend.cex = legend_font_size
      )

      image(
        t(d2prob[genes_select, ]),
        main = paste(celltypes[j], "regulators"),
        col = custom_colors,
        axes = FALSE
      )
      axis(
        1,
        at = seq(0, 1, length.out = ncol(d2prob)),
        labels = colnames(d2prob),
        cex.axis = font_size
      )
      axis(
        2,
        at = seq(0, 1, length.out = ntop),
        labels = genes_select,
        las = 2,
        cex.axis = font_size
      )

      fields::image.plot(
        add = TRUE,
        legend.only = TRUE,
        zlim = range(d2prob),
        col = custom_colors,
        legend.width = legend_width,
        legend.cex = legend_font_size
      )

      plot(
        mydist[genes_select],
        seq_along(genes_select),
        type = "o",
        main = paste("Regulator rewiring score:
                     ", celltypes[i], "vs", celltypes[j]),
        col = line_color,
        las = 1,
        xlab = NA,
        ylab = NA,
        cex = font_size,
        yaxt = "n"
      )
      axis(
        2,
        at = seq_along(genes_select),
        labels = genes_select,
        las = 2,
        cex.axis = font_size
      )

      barplot(
        data1_reg_degree[genes_select],
        horiz = TRUE,
        main = paste(celltypes[i], "regulators"),
        xlim = c(0, max(data1_reg_degree[genes_select]) + 0.5),
        col = bar_color,
        las = 1,
        cex.names = font_size,
        names.arg = genes_select
      )
      barplot(
        data2_reg_degree[genes_select],
        horiz = TRUE,
        main = paste(celltypes[j], "regulators"),
        xlim = c(0, max(data2_reg_degree[genes_select]) + 0.5),
        col = bar_color,
        las = 1,
        cex.names = font_size,
        names.arg = genes_select
      )
      dev.off()
    }
  }

  return(genes_select_tables[, -1])
}
