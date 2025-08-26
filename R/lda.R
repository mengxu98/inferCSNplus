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

  networks_matrix_list <- purrr::map(
    network_list,
    .f = function(x) {
      colnames(x) <- c("row", "col", "value")
      network_matrix <- thisutils::table_to_matrix(x)
      network_matrix <- abs(network_matrix)
      rowsum <- as.numeric(rowSums(network_matrix))
      idrms <- as.numeric(which(rowsum == 0))
      if (length(idrms) != 0) {
        network_matrix <- network_matrix[-idrms, ]
      }
      as.data.frame(network_matrix)
    }
  )

  tfs_list <- purrr::map(
    networks_matrix_list,
    .f = function(x) {
      regulators <- rownames(x)
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

  if (combine_method == "intersect") {
    networks_matrix <- purrr::map_dfr(
      networks_matrix_list,
      .f = function(x) {
        x[tfs_list_intersect, ]
      }
    )
  } else if (combine_method == "c") {
    networks_matrix <- purrr::map_dfr(
      networks_matrix_list,
      .f = function(x) {
        x
      }
    )
  }
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

  networks_matrix <- abs(networks_matrix)
  rowsum <- rowSums(networks_matrix)
  idrms <- which(rowsum == 0)

  set.seed(seed)
  if (is.null(path)) {
    path <- "lda/"
  }
  outdir <- file.path(path, paste0("k", k))
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (binary) {
    networks_matrix[networks_matrix >= binary_threshold] <- 1
    networks_matrix[networks_matrix < binary_threshold] <- 0
    lda_result <- topicmodels::LDA(networks_matrix, k)
  } else {
    # lda_result <- topicmodels::LDA(ceiling(networks_matrix * 100), k)
    lda_result <- topicmodels::LDA(ceiling(networks_matrix), k)
  }

  save(
    lda_result,
    file = file.path(outdir, paste0("lda_model_k", k, ".RData"))
  )

  # Plotting Document-Topic and Topic-Word Probabilities
  doc_prob <- modeltools::posterior(lda_result)$topics
  rownames(doc_prob) <- regulators
  word_prob <- modeltools::posterior(lda_result)$terms

  grDevices::pdf(
    file = file.path(outdir, paste0("lda_k", k, "_allcells", ".pdf")),
    width = 8,
    height = 5
  )
  graphics::par(mfrow = c(1, 2))
  image(
    t(doc_prob),
    main = "Document-Topic Probabilities",
    xlab = "Topics",
    ylab = "TFs * cells",
    col = custom_colors,
    axes = FALSE
  )
  graphics::axis(
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
  graphics::axis(
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
  grDevices::dev.off()

  docprob <- doc_prob
  topicid <- apply(docprob, 1, which.max)
  tforder <- order(topicid)

  grDevices::pdf(
    file = file.path(outdir, paste0("lda_k", k, "_tf-topicPercell.pdf")),
    width = 15,
    height = 10
  )
  graphics::par(mfrow = c(2, 4))
  for (i in seq_along(celltypes)) {
    id <- start_index[i]:end_index[i]
    idrm <- intersect(idrms, id)
    idkp <- setdiff(id, idrm)
    idkp <- intersect(idkp, seq_len(nrow(docprob)))

    d1prob <- docprob[idkp, ]
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
    graphics::axis(
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
  grDevices::dev.off()

  # Save genes and topic ids
  TopicWordProb <- apply(t(word_prob), 1, max)
  geneid <- apply(t(word_prob), 1, which.max)
  utils::write.table(
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
    utils::write.table(
      data.frame(names(TopicProb), TFtopicid, TopicProb),
      file = tf_topic_file,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE
    )

    d1prob <- docprob[idkp, ]
    tfid <- apply(d1prob, 1, which.max)
    tforder <- order(tfid)

    grDevices::pdf(
      file = file.path(outdir, paste0("lda_k", k, "_tf-topicPercell.pdf")),
      width = 15,
      height = 10
    )
    graphics::par(mfrow = c(2, 4))
    for (i in seq_along(celltypes)) {
      id <- start_index[i]:end_index[i]
      idrm <- intersect(idrms, id)
      idkp <- setdiff(id, idrm)

      d1prob <- docprob[idkp, ]
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
      graphics::axis(
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
    grDevices::dev.off()
  }

  genes_select_tables <- as.data.frame(
    matrix(data = 1, nrow = ntop)
  )
  # Comparing TF-topic distribution between cell types
  for (i in seq_along(celltypes)) {
    id1 <- start_index[i]:end_index[i]
    d1prob <- docprob[id1, ]
    data1_reg <- networks_matrix[id1, ]
    rownames(data1_reg) <- regulators[id1]
    data1_reg_degree <- rowSums(data1_reg)

    # if (i < length(celltypes)) {
    #   j <- i + 1
    # } else {
    #   break
    # }

    for (j in (i + 1):length(celltypes)) {
      if (i == length(celltypes)) {
        j <- 1
      }
      id2 <- start_index[j]:end_index[j]
      d2prob <- docprob[id2, ]
      data2_reg <- networks_matrix[id2, ]
      rownames(data2_reg) <- regulators[id2]
      data2_reg_degree <- rowSums(data2_reg)

      divergence_res <- (
        kl_divergence(d1prob, d2prob) + kl_divergence(d2prob, d1prob)
      ) / 2

      divergence_res <- divergence_res[order(divergence_res, decreasing = TRUE)]

      genes_select <- names(divergence_res)[1:ntop]

      genes_select_table <- data.frame(
        genes_select
      )
      names(genes_select_table) <- paste0(celltypes[i], "_", celltypes[j])
      genes_select_tables <- cbind(genes_select_tables, genes_select_table)

      grDevices::pdf(
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

      graphics::par(mfrow = c(1, 5))

      image(
        t(d1prob[genes_select, ]),
        main = paste(celltypes[i], "regulators"),
        col = custom_colors,
        axes = FALSE
      )
      graphics::axis(
        1,
        at = seq(0, 1, length.out = ncol(d1prob)),
        labels = colnames(d1prob),
        cex.axis = font_size
      )
      graphics::axis(
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
      graphics::axis(
        1,
        at = seq(0, 1, length.out = ncol(d2prob)),
        labels = colnames(d2prob),
        cex.axis = font_size
      )
      graphics::axis(
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
        divergence_res[genes_select],
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
      graphics::axis(
        2,
        at = seq_along(genes_select),
        labels = genes_select,
        las = 2,
        cex.axis = font_size
      )

      graphics::barplot(
        data1_reg_degree[genes_select],
        horiz = TRUE,
        main = paste(celltypes[i], "regulators"),
        xlim = c(0, max(data1_reg_degree[genes_select]) + 0.5),
        col = bar_color,
        las = 1,
        cex.names = font_size,
        names.arg = genes_select
      )
      graphics::barplot(
        data2_reg_degree[genes_select],
        horiz = TRUE,
        main = paste(celltypes[j], "regulators"),
        xlim = c(0, max(data2_reg_degree[genes_select]) + 0.5),
        col = bar_color,
        las = 1,
        cex.names = font_size,
        names.arg = genes_select
      )
      grDevices::dev.off()
    }
  }

  return(genes_select_tables[, -1])
}
