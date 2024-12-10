#' @title rank umap
#'
#' @param umap input
#'
#' @return output
#' @export
vector_rank_pca <- function(umap) {
  pca_out <- gmodels::fast.prcomp(
    apply(umap, 2, rank),
    retx = TRUE,
    center = TRUE,
    scale. = TRUE,
    tol = NULL
  )

  return(pca_out$x)
}

#' Title
#'
#' @param vector vector
#' @param N N
#' @param plot plot
#' @param colors colors
#'
#' @return output
#' @export
vector_build_grid <- function(
    vector,
    N = 30,
    plot = TRUE,
    colors = "grey70") {
  vector_e <- vector
  delta <- 0.000001

  if (plot == TRUE) {
    graphics::plot(
      vector_e,
      col = colors,
      pch = 16,
      cex = 0.2,
      xlab = "",
      ylab = "",
      axes = FALSE
    )
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black"
    )
  }

  X <- vector[, 1]
  Y <- vector[, 2]

  this_step_x <- (max(X) - min(X) - delta) / N
  this_step_y <- (max(Y) - min(Y) - delta) / N
  NUM_CUT <- 1

  index_list <- list()
  center_list <- list()

  this_x <- min(X)
  while (this_x < max(X)) {
    this_y <- min(Y)
    while (this_y < max(Y)) {
      this_in_index <- which(
        vector_e[, 1] >= this_x & vector_e[, 1] < this_x + this_step_x &
          vector_e[, 2] >= this_y & vector_e[, 2] < this_y + this_step_y
      )
      this_center <- c(this_x + this_step_x / 2, this_y + this_step_y / 2)
      if (length(this_in_index) >= NUM_CUT) {
        index_list <- c(index_list, list(this_in_index))
        center_list <- c(center_list, list(this_center))
        if (plot == TRUE) {
          graphics::points(
            this_center[1],
            this_center[2],
            col = "black",
            pch = 16,
            cex = 0.5
          )
        }
      }
      this_y <- this_y + this_step_y
    }
    this_x <- this_x + this_step_x
  }

  out <- list()
  out$vector <- vector_e
  out$index_list <- index_list
  out$center_list <- center_list
  out$this_step_x <- this_step_x
  out$this_step_y <- this_step_y

  return(out)
}

#' Title
#'
#' @param out out
#' @param CUT CUT
#' @param plot plot
#' @param colors colors
#'
#' @return output
#' @export
vector_build_net <- function(
    out,
    CUT = 1,
    plot = TRUE,
    colors = "grey70") {
  vector <- out$vector
  index_list <- out$index_list
  center_list <- out$center_list
  this_step_x <- out$this_step_x
  this_step_y <- out$this_step_y
  delta <- 0.00001

  if (plot == TRUE) {
    graphics::plot(
      vector,
      col = colors,
      pch = 16,
      cex = 0.2,
      xlab = "",
      ylab = "",
      axes = FALSE
    )
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black"
    )
  }

  center_vector <- c()
  i <- 1
  while (i <= length(center_list)) {
    center_vector <- cbind(center_vector, center_list[[i]])
    i <- i + 1
  }
  center_vector <- t(center_vector)
  out$center_vector <- center_vector

  p1 <- c()
  p2 <- c()

  CNUM <- length(center_list)
  i <- 1
  while (i <= CNUM) {
    this_p1_loc <- center_list[[i]]
    this_p1 <- paste0("P", as.character(i))

    used_j <- which(
      (abs(center_vector[, 1] - this_p1_loc[1]) <= this_step_x + delta) &
        (abs(center_vector[, 2] - this_p1_loc[2]) <= this_step_y + delta)
    )

    for (j in used_j) {
      this_p2 <- paste0("P", as.character(j))

      this_p2_loc <- center_list[[j]]

      if (length(index_list[[i]]) >= CUT &&
        length(index_list[[j]]) >= CUT &&
        this_p1 != this_p2) {
        p1 <- c(p1, this_p1)
        p2 <- c(p2, this_p2)
        if (plot == TRUE) {
          graphics::segments(
            x0 = this_p1_loc[1], x1 = this_p2_loc[1],
            y0 = this_p1_loc[2], y1 = this_p2_loc[2],
            col = "black"
          )
        }
      }
    }
    i <- i + 1
  }

  tag <- c()
  i <- 1
  while (i <= length(p1)) {
    this_p1 <- p1[i]
    this_p2 <- p2[i]
    sorted_pair <- sort(c(this_p1, this_p2))
    this_tag <- paste0(sorted_pair[1], "|", sorted_pair[2])
    tag <- c(tag, this_tag)
    i <- i + 1
  }
  tag <- unique(tag)
  p1 <- c()
  p2 <- c()
  i <- 1
  while (i <= length(tag)) {
    this_p1 <- strsplit(tag[i], "\\|")[[1]][1]
    this_p2 <- strsplit(tag[i], "\\|")[[1]][2]
    p1 <- c(p1, this_p1)
    p2 <- c(p2, this_p2)
    i <- i + 1
  }

  out$p1 <- p1
  out$p2 <- p2
  NET <- cbind(p1, p2)
  g <- igraph::make_graph(t(NET), directed = FALSE)
  ALLNAME <- paste0("P", 1:CNUM)
  ADD <- ALLNAME[which(!ALLNAME %in% igraph::as_ids(igraph::V(g)))]
  g <- g + ADD

  dist <- igraph::distances(
    g,
    v = igraph::V(g),
    to = igraph::V(g),
    mode = c("all")
  )

  DIST.NUM <- as.numeric(
    stringr::str_replace(colnames(dist), "P", "")
  )
  dist <- dist[order(DIST.NUM), order(DIST.NUM)]

  CPT <- igraph::components(g)
  MAXC <- which(CPT$csize == max(CPT$csize))[1]

  used_name <- names(which(CPT$membership == MAXC))
  used <- as.numeric(
    stringr::str_replace(used_name, "P", "")
  )

  used_name <- used_name[order(used)]
  used <- used[order(used)]
  if (plot == TRUE) {
    graphics::points(out$center_vector[used, ], col = "red", pch = 16, cex = 0.5)
  }
  USED_INDEX <- c()
  i <- 1
  while (i <= length(used)) {
    USED_INDEX <- c(USED_INDEX, index_list[[used[i]]])
    i <- i + 1
  }

  out$graph <- g
  out$dist <- dist
  out$used <- used
  out$used_name <- used_name
  out$USED_INDEX <- USED_INDEX

  return(out)
}

#' Title
#'
#' @param out out
#' @param umap umap
#' @param plot plot
#'
#' @return output
#' @export
vector_build_value <- function(out, umap, plot = TRUE) {
  VALUE.OUT <- vector.calValue(umap)

  out$value <- VALUE.OUT$value
  out$umap <- umap
  out$PCA.RC <- VALUE.OUT$PCA.RC

  if (plot == TRUE) {
    vector.showValue(out)
  }

  return(out)
}

#' Title
#'
#' @param out out
#' @param plot plot
#'
#' @return output
#' @export
vector_grid_value <- function(out, plot = TRUE) {
  index_list <- out$index_list
  value <- out$value
  plot <- plot
  used <- out$used

  center_value <- c()
  i <- 1
  while (i <= length(index_list)) {
    this_value <- mean(value[index_list[[i]]])
    center_value <- c(center_value, this_value)
    i <- i + 1
  }
  # center_vector <- out$center_vector

  value_norm <- normalization(value)
  VALUE.COL <- vector_vcol(value_norm, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  value <- center_value
  value_norm <- normalization(value)
  colors <- vector_vcol(value_norm, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  if (plot == TRUE) {
    graphics::plot(
      out$vector,
      col = "grey80",
      pch = 16,
      cex = 0.5,
      xlab = "",
      ylab = "",
      axes = FALSE
    )
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black"
    )
    # graphics::points(out$center_vector,col=colors,pch=16)
    graphics::points(
      out$center_vector[used, ],
      col = colors[used],
      pch = 15,
      cex = 1.5
    )
  }
  out$center_value <- center_value
  out$ORIG.CENTER.COL <- colors
  out$VALUE.COL <- VALUE.COL

  return(out)
}

#' Title
#'
#' @param out out
#' @param up_value up_value
#' @param plot plot
#'
#' @return output
#' @export
vector_auto_center <- function(
    out,
    up_value = 0.9,
    plot = TRUE) {
  dist <- out$dist
  used <- out$used
  used_name <- out$used_name
  center_value <- out$center_value
  center_vector <- out$center_vector
  # center_index <- out$center_index
  index_list <- out$index_list

  used_center_value <- center_value[used]
  # USED_CENTER_VEC <- center_vector[used, ]

  high <- used[which(used_center_value >= stats::quantile(used_center_value, up_value))]
  high_name <- used_name[which(used_center_value >= stats::quantile(used_center_value, up_value))]
  # LOW_NAME <- used_name[which(used_center_value < stats::quantile(used_center_value, up_value))]
  # plot(out$center_vector[used,])
  # graphics::points(center_vector[high,], col='red',pch=16)
  # plot(center_vector[high,], col='red')

  sub_graph <- igraph::induced_subgraph(out$graph, high_name)
  sub_cpt <- igraph::components(sub_graph)
  cluster <- list()
  length <- c()
  pch <- rep(1, length(high))
  dist_cor <- c()
  dist_mean <- c()

  i <- 1
  while (i <= sub_cpt$no) {
    this_name <- names(which(sub_cpt$membership == i))
    this_index <- as.numeric(stringr::str_replace(this_name, "P", ""))
    pch[which(high %in% this_index)] <- as.character(i)
    length <- c(length, length(this_index))
    if (length(this_index) == 1) {
      this_dist <- dist[used, this_index]
    } else {
      this_dist <- apply(dist[used, this_index], 1, mean)
    }
    this_cor <- stats::cor(this_dist, center_value[used], method = "spearman")
    dist_cor <- c(dist_cor, this_cor)

    dist_mean <- c(dist_mean, mean(this_dist))
    cluster <- c(cluster, list(this_index))
    i <- i + 1
  }

  # select=which(dist_cor==min(dist_cor))[1]
  # select=which(rank(dist_mean)*rank(dist_cor) == min(rank(dist_mean)*rank(dist_cor)) )

  tmp <- 10^length - dist_cor
  # tmp=rank(-dist_cor) * rank(length )
  select <- which(tmp == max(tmp))[1]
  # select=which(length==max(length))[1]
  # select=which( rank(-dist_cor) * rank(length) == max( rank(-dist_cor) * rank(length)  ) )[1]

  summit <- cluster[[select]]
  # graphics::points(out$center_vector[cluster[[9]],],pch=16,col='red')

  score <- c()
  PS <- c()
  i <- 1
  while (i <= length(used_name)) {
    this_name <- used_name[i]
    this_dist <- dist[
      which(colnames(dist) == this_name),
      which(rownames(dist) %in% paste0("P", summit))
    ]
    this_score <- min(this_dist)
    score <- c(score, this_score)
    PS <- c(PS, this_score)
    i <- i + 1
  }
  score <- max(score) - score

  value <- score
  # plot(out$vector, col='grey70',pch=16)
  value_norm <- normalization(value)
  colors <- vector_vcol(
    value_norm,
    c(0, 0.5, 1),
    c("#009FFF", "#FFF200", "#ec2F4B")
  ) # too strong
  # colors <- vector_vcol(value_norm, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ee9ca7"))

  out$colors <- rep("grey70", nrow(out$vector))
  out$ORIG.COL <- rep("grey70", nrow(out$vector))
  out$p_score <- rep(0, nrow(out$vector))
  out$P.PS <- rep(NA, nrow(out$vector))
  i <- 1
  while (i <= length(used)) {
    this_index <- index_list[[used[i]]]
    out$colors[this_index] <- colors[i]
    out$ORIG.COL[this_index] <- out$ORIG.CENTER.COL[used][i]
    out$p_score[this_index] <- score[i]
    out$P.PS[this_index] <- PS[i]
    i <- i + 1
  }

  if (plot == TRUE) {
    graphics::plot(
      out$vector,
      col = out$ORIG.COL,
      pch = 16,
      cex = 0.5,
      xlab = "",
      ylab = "",
      axes = FALSE
    )
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black"
    )
    graphics::text(
      center_vector[high, 1],
      center_vector[high, 2],
      labels = pch,
      cex = 1,
      pos = 2
    )
    graphics::points(
      center_vector[high, 1],
      center_vector[high, 2],
      col = "black",
      pch = 16,
      cex = 1
    )
    graphics::points(
      center_vector[summit, 1],
      center_vector[summit, 2],
      col = "black",
      pch = 16,
      cex = 1.5
    )
    graphics::points(
      center_vector[summit, 1],
      center_vector[summit, 2],
      col = "red",
      pch = 16,
      cex = 1
    )
  }

  out$score <- score
  out$summit <- summit
  out$cluster <- cluster
  out$length <- length
  out$pch <- pch
  out$dist_cor <- dist_cor
  out$PS <- PS
  # out$dist_mean=dist_mean

  return(out)
}

#' Title
#'
#' @param out out
#' @param P P
#' @param plot plot
#' @param colors colors
#' @param OL OL
#' @param arrow_length2 arrow_length2
#' @param cex_value cex_value
#' @param arrow_width arrow_width
#' @param BD BD
#' @param arrow_color arrow_color
#' @param plot.SUMMIT plot.SUMMIT
#'
#' @return output
#' @export
vector_draw_arrow <- function(
    out,
    P = 0.9,
    plot = TRUE,
    colors = "grey70",
    OL = 1.5,
    arrow_length2 = 70,
    cex_value = 0.5,
    arrow_width = 1,
    BD = TRUE,
    arrow_color = "grey30",
    plot.SUMMIT = TRUE) {
  used <- out$used
  dist <- out$dist
  ALL_VEC <- out$vector
  used_name <- out$used_name
  USED_CENTER_VEC <- out$center_vector[used, ]
  USED_DIST <- out$dist[which(rownames(dist) %in% used_name), which(rownames(dist) %in% used_name)]
  score <- out$score

  one <- min(stats::dist(USED_CENTER_VEC)) * OL

  DIV <- 1 / P

  if (plot == TRUE) {
    if (BD == TRUE) {
      graphics::plot(
        ALL_VEC,
        col = colors,
        pch = 16,
        cex = cex_value,
        xlim = c(min(ALL_VEC[, 1]) - one, max(ALL_VEC[, 1]) + one),
        ylim = c(min(ALL_VEC[, 2]) - one, max(ALL_VEC[, 2]) + one),
        xlab = "",
        ylab = "",
        axes = FALSE
      )
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NA,
        border = "black"
      )
    } else {
      graphics::plot(
        ALL_VEC,
        col = colors,
        pch = 16,
        cex = cex_value,
        xlim = c(min(ALL_VEC[, 1]) - one, max(ALL_VEC[, 1]) + one),
        ylim = c(min(ALL_VEC[, 2]) - one, max(ALL_VEC[, 2]) + one),
        yaxt = "n",
        axes = F,
        xlab = "",
        ylab = "",
        axes = FALSE
      )
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NA,
        border = "black"
      )
    }
  }
  N.SCORE <- normalization(score)
  SCORE.COL <- vector_vcol(
    N.SCORE,
    c(0, 0.5, 1),
    c("#009FFF", "#FFF200", "#ec2F4B")
  )

  A1_VEC <- c()
  A2_VEC <- c()
  A_LENGTH <- c()

  i <- 1
  while (i <= length(used)) {
    this_p1_loc <- USED_CENTER_VEC[i, ]

    vector_list <- cbind(USED_CENTER_VEC[, 1] - this_p1_loc[1], USED_CENTER_VEC[, 2] - this_p1_loc[2])
    vector_list_norm <- t(apply(vector_list, 1, .norm_one, one))

    vector_weight_1 <- DIV^-(rank(USED_DIST[i, ]) - 1) # equals to P ^ (rank(USED_DIST[i,])-1)
    vector_weight_2 <- score[i] - score # PS - PS[i]; "PS_i - PS_j" in  the paper; In paper, PS[i] is "PS_j".

    vector_weight <- vector_weight_1 * vector_weight_2
    vector_weight <- vector_weight / sum(abs(vector_weight))

    final_vec <- t(vector_list_norm) %*% vector_weight

    this_p2_loc <- c(this_p1_loc[1] + final_vec[1], this_p1_loc[2] + final_vec[2])

    arrow_length <- grDevices::dev.size()[1] / arrow_length2 * sqrt(sum(final_vec^2)) / one # * sqrt(sum(final_vec^2)) #0.25
    if (plot == TRUE) {
      graphics::arrows(
        x0 = this_p1_loc[1],
        y0 = this_p1_loc[2],
        x1 = this_p2_loc[1],
        y1 = this_p2_loc[2],
        lwd = arrow_width,
        length = arrow_length,
        col = arrow_color
      )
    }
    A1_VEC <- cbind(A1_VEC, this_p1_loc)
    A2_VEC <- cbind(A2_VEC, this_p2_loc)
    A_LENGTH <- c(A_LENGTH, arrow_length)
    i <- i + 1
  }

  if (plot == TRUE & plot.SUMMIT == TRUE) {
    X1 <- min(out$center_vector[out$summit, 1]) - one / 10
    X2 <- max(out$center_vector[out$summit, 1]) + one / 10
    Y1 <- min(out$center_vector[out$summit, 2]) - one / 10
    Y2 <- max(out$center_vector[out$summit, 2]) + one / 10

    graphics::rect(
      xleft = X1,
      ybottom = Y1,
      xright = X2,
      ytop = Y2,
      angle = 45,
      col = NA,
      border = "#009FFF",
      lty = 1,
      lwd = 2
    )
  }
  A1_VEC <- t(A1_VEC)
  A2_VEC <- t(A2_VEC)
  out$A1_VEC <- A1_VEC
  out$A2_VEC <- A2_VEC
  out$A_LENGTH <- A_LENGTH
  out$A_COL <- SCORE.COL

  return(out)
}

# Other functions
select_Seurat <- function(
    object,
    type = "select") {
  type <- match.arg(type, c("select", "unselect"))
  p <- Seurat::DimPlot(
    object,
    reduction = "umap",
    pt.size = 0.5
  )
  selected_cells <- Seurat::CellSelector(plot = p)
  if (type == "select") {
    return(selected_cells)
  } else if (type == "unselect") {
    return(setdiff(colnames(object), selected_cells))
  } else {
    return()
  }
}

vector.AddMetaByCell.Seurat <- function(
    object,
    used_cells) {
  select <- rep("NO", ncol(object))
  select[which(colnames(object) %in% used_cells)] <- "YES"
  object@meta.data$select <- select
  return(object)
}

vector.pca.Seurat <- function(
    object,
    method = "prcomp",
    random = FALSE,
    RN = 1000,
    CUT = 0.5,
    seed = 1) {
  set.seed(seed = seed)

  D <- as.matrix(object@assays$RNA@scale.data)
  if (random) {
    R_INDEX <- sample(1:ncol(D), RN, replace = TRUE)
    D_raw <- D
    D <- D[, R_INDEX]
  }

  pca_out <- switch(
    EXPR = method,
    "fast.prcomp" = gmodels::fast.prcomp(
      t(D),
      retx = TRUE,
      center = FALSE,
      scale. = FALSE,
      tol = NULL
    ),
    "prcomp" = stats::prcomp(t(D))
  )

  EXP <- (cumsum(pca_out$sdev^2) / sum(pca_out$sdev^2))
  N <- min(which(cumsum(pca_out$sdev^2) / sum(pca_out$sdev^2) > CUT))

  if (random) {
    pred_pca <- t(D_raw) %*% pca_out$rotation
  }

  pca_out$EXP <- EXP
  pca_out$CUT <- CUT
  pca_out$N <- N
  if (random) {
    pca_out$RN <- RN
    pca_out$R_INDEX <- R_INDEX
    pca_out$pred_pca <- pred_pca
  }

  return(pca_out)
}

vector.lcol <- function(tag) {
  tag <- as.factor(tag)
  my_color_palette <- scales::hue_pal()(length(unique(tag)))
  colors <- my_color_palette[tag]
  return(colors)
}

vector_vcol <- function(value, CV, CN) {
  CRF <- circlize::colorRamp2(CV, CN)
  colors <- CRF(value)
  return(colors)
}

vector.showValue <- function(out) {
  vector <- out$vector
  value <- out$value
  value_norm <- normalization(value)
  colors <- vector_vcol(
    value_norm,
    c(0, 0.5, 1),
    c("#009FFF", "#FFF200", "#ec2F4B")
  )
  graphics::plot(
    vector,
    col = colors,
    pch = 16,
    cex = 0.5,
    xlab = "",
    ylab = "",
    axes = FALSE
  )
  graphics::rect(
    graphics::par("usr")[1],
    graphics::par("usr")[3],
    graphics::par("usr")[2],
    graphics::par("usr")[4],
    col = NA,
    border = "black"
  )

  return(out)
}

vector.calValue <- function(umap) {
  out <- list()
  PCA.RC <- apply(apply(umap, 2, rank), 2, normalization)
  PCA.RC <- abs(PCA.RC - 0.5)
  value <- apply(PCA.RC, 1, mean)
  out$value <- value
  out$PCA.RC <- PCA.RC
  return(out)
}

vector.nonCenter <- function(out) {
  out$score <- out$center_value[out$used]
  out$colors <- out$VALUE.COL
  return(out)
}

.norm_one <- function(x, one = 1) {
  if (stats::var(x) != 0) {
    x <- x / sqrt(sum(x^2)) * one
  }
  return(x)
}

# Manually do something
vector.selectPoint <- function(vector, cex_value = 0.5) {
  points <- gatepoints::fhs(vector, pch = 16, col = "red3", cex = cex_value, mark = TRUE)
  return(points)
}

vector.selectCenter <- function(out, plot = TRUE) {
  dist <- out$dist
  used <- out$used
  used_name <- out$used_name
  center_vector <- out$center_vector
  index_list <- out$index_list

  graphics::plot(out$vector, col = "grey80", pch = 16, cex = 0.5)
  graphics::points(out$center_vector[used, ], col = out$ORIG.CENTER.COL[used], pch = 16, cex = 1)
  SELECT_NAME <- vector.selectPoint(out$center_vector[used, ], cex_value = 1)
  summit <- used[as.numeric(SELECT_NAME)]

  score <- c()
  i <- 1
  while (i <= length(used_name)) {
    this_name <- used_name[i]
    this_dist <- dist[which(colnames(dist) == this_name), which(rownames(dist) %in% paste0("P", summit))]
    this_dist <- this_dist
    # this_value=center_value[used]
    this_score <- min(this_dist) # sum(rank(-this_dist) * rank(this_value))
    # this_cor=cor(-this_value, this_dist)#,method='spearman')
    score <- c(score, this_score)
    i <- i + 1
  }

  PS <- score

  score <- max(score) - score

  value <- score
  value_norm <- normalization(value)
  colors <- vector_vcol(value_norm, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  out$colors <- rep("grey70", nrow(out$vector))
  out$ORIG.COL <- rep("grey70", nrow(out$vector))
  out$p_score <- rep(0, nrow(out$vector))
  out$P.PS <- rep(NA, nrow(out$vector))
  i <- 1
  while (i <= length(used)) {
    this_index <- index_list[[used[i]]]
    out$colors[this_index] <- colors[i]
    out$ORIG.COL[this_index] <- out$ORIG.CENTER.COL[used][i]
    out$p_score[this_index] <- score[i]
    out$P.PS[this_index] <- PS[i]
    i <- i + 1
  }

  if (plot == TRUE) {
    graphics::plot(out$vector, col = out$colors, pch = 16, cex = 0.5)
    graphics::points(center_vector[summit, ], col = "black", pch = 16, cex = 1.5)
    graphics::points(center_vector[summit, ], col = "red", pch = 16, cex = 1)
  }

  out$score <- score
  out$PS <- PS
  out$summit <- summit

  return(out)
}

vector.reDrawArrow <- function(out, colors = "grey70") {
  A1_VEC <- out$A1_VEC
  A2_VEC <- out$A2_VEC
  A_LENGTH <- out$A_LENGTH
  vector <- out$vector

  graphics::plot(x = vector[, 1], y = vector[, 2], col = colors, cex = 0.5, pch = 16)
  i <- 1
  while (i <= length(A_LENGTH)) {
    graphics::arrows(
      x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
      x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
      lwd = 2, length = A_LENGTH[i],
      col = "black"
    )
    i <- i + 1
  }

  return(out)
}


vector.selectRegion <- function(out) {
  p_score <- out$p_score
  A1_VEC <- out$A1_VEC
  A2_VEC <- out$A2_VEC
  A_LENGTH <- out$A_LENGTH
  vector <- out$vector
  P.PS <- out$P.PS

  index_list <- out$index_list
  used <- out$used
  SELECT_NAME <- vector.selectPoint(vector, cex_value = 0.1)
  SELECT_INDEX <- which(rownames(vector) %in% SELECT_NAME)
  SELECT_INDEX <- SELECT_INDEX[which(SELECT_INDEX %in% out$USED_INDEX)]

  A_USED <- c()
  i <- 1
  while (i <= length(used)) {
    if (length(which(index_list[[used[i]]] %in% SELECT_INDEX)) > 0) {
      A_USED <- c(A_USED, i)
    }

    i <- i + 1
  }

  # Draw new
  colors <- out$colors
  colors[which(!rownames(vector) %in% SELECT_NAME)] <- "grey70"
  graphics::plot(x = vector[, 1], y = vector[, 2], col = colors, cex = 0.5, pch = 16)
  i <- 1
  while (i <= length(A_LENGTH)) {
    graphics::arrows(
      x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
      x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
      lwd = 2, length = A_LENGTH[i],
      col = "black"
    )
    i <- i + 1
  }

  for (i in A_USED) {
    graphics::arrows(
      x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
      x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
      lwd = 2, length = A_LENGTH[i],
      col = "red"
    )
  }

  out$A_USED <- A_USED
  out$SELECT_NAME <- SELECT_NAME
  out$SELECT_INDEX <- SELECT_INDEX
  out$SELECT_SCORE <- p_score[SELECT_INDEX]
  out$SELECT_PS <- P.PS[SELECT_INDEX]

  return(out)
}

# vector.reDrawRegion <- function(out) {
#   p_score <- out$p_score
#   A1_VEC <- out$A1_VEC
#   A2_VEC <- out$A2_VEC
#   A_LENGTH <- out$A_LENGTH
#   vector <- out$vector
#   P.PS <- out$P.PS
#   index_list <- out$index_list
#   used <- out$used
#   SELECT_NAME <- out$SELECT_NAME # vector.selectPoint(vector,cex_value=0.1)
#   SELECT_INDEX <- which(rownames(vector) %in% SELECT_NAME)

#   A_USED <- c()
#   i <- 1
#   while (i <= length(used)) {
#     if (length(which(index_list[[used[i]]] %in% SELECT_INDEX)) > 0) {
#       A_USED <- c(A_USED, i)
#     }

#     i <- i + 1
#   }

#   # Draw new
#   colors <- out$colors
#   colors[which(!rownames(vector) %in% SELECT_NAME)] <- "grey70"
#   graphics::plot(x = vector[, 1], y = vector[, 2], col = colors, cex = 0.5, pch = 16)

#   i <- 1
#   while (i <= length(A_LENGTH)) {
#     graphics::arrows(
#       x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
#       x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
#       lwd = 2, length = A_LENGTH[i],
#       col = "black"
#     )
#     i <- i + 1
#   }

#   for (i in A_USED) {
#     graphics::arrows(
#       x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
#       x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
#       lwd = 2, length = A_LENGTH[i],
#       col = "red"
#     )
#   }
#   out$A_USED <- A_USED
#   out$SELECT_NAME <- SELECT_NAME
#   out$SELECT_INDEX <- SELECT_INDEX
#   out$SELECT_SCORE <- p_score[SELECT_INDEX]
#   out$SELECT_PS <- P.PS[SELECT_INDEX]

#   return(out)
# }

# vector.RPPCA <- function(umap) {
#   R.PCA <- apply(umap, 2, rank)
#   D <- t(R.PCA)
#   pca_out <- gmodels::fast.prcomp(
#     t(D),
#     retx = TRUE,
#     center = FALSE,
#     scale. = FALSE,
#     tol = NULL
#   )
#   PPCA <- pca_out$x
#   return(PPCA)
# }


# vector.medCurv <- function(
#     umap,
#     MAX = 1000) {
#   MS <- vector.calValue(umap)$PCA.RC
#   COR <- cor(MS)
#   MED <- c()
#   ALL <- c()
#   i <- 1
#   while (i <= ncol(COR) && i <= MAX) {
#     this_col <- COR[1:i, i]
#     ALL <- c(ALL, this_col)
#     ALL[which(ALL == 1)] <- NA
#     this_med <- median(ALL, na.rm = TRUE)
#     MED <- c(MED, this_med)
#     if (i %% 100 == 1) {
#       print(i)
#     }
#     i <- i + 1
#   }
#   MED[1] <- MED[2]
#   return(MED)
# }

# 2020.01.03

# vector.smoothOut <- function(X, Z) {
#   Y <- X
#   Y[order(Z)] <- X[order(Z)] - stats::smooth(X[order(Z)])
#   return(Y)
# }

# vector.regressOut <- function(X, Z) {
#   FIT <- stats::lm(X ~ Z)
#   Y <- X - predict(FIT)
#   return(Y)
# }

# vector.removeOut <- function(X) {
#   Q1 <- stats::quantile(X, 0.25)
#   Q3 <- stats::quantile(X, 0.75)
#   IQR <- Q3 - Q1
#   LW <- Q1 - 1.5 * IQR
#   up_value <- Q3 + 1.5 * IQR
#   X[which(X > up_value)] <- up_value
#   X[which(X < LW)] <- LW

#   return(X)
# }

# 2020.1.12
# .normNew <- function(x) {
#   pos_index <- which(x > 0)
#   neg_index <- which(x < 0)
#   x[pos_index] <- rank(x[pos_index])
#   x[neg_index] <- rank(-x[neg_index])

#   return(x)
# }

# vector.calValueNew <- function(umap) {
#   out <- list()
#   PCA.RC <- apply(umap, 2, .normNew)
#   value <- apply(PCA.RC, 1, mean)
#   out$value <- value
#   out$PCA.RC <- PCA.RC

#   return(out)
# }

# vector.getValueNew <- function(out, umap, plot = TRUE) {
#   VALUE.OUT <- vector.calValueNew(umap)
#   out$value <- VALUE.OUT$value
#   out$umap <- umap
#   out$PCA.RC <- VALUE.OUT$PCA.RC
#   if (plot == TRUE) {
#     vector.showValue(out)
#   }

#   return(out)
# }


# 20201030
# vector.autoCenterCor <- function(
#     out,
#     up_value = 0.9,
#     plot = TRUE) {
#   dist <- out$dist
#   used <- out$used
#   used_name <- out$used_name
#   center_value <- out$center_value
#   center_vector <- out$center_vector
#   index_list <- out$index_list

#   used_center_value <- center_value[used]

#   high <- used[which(used_center_value >= stats::quantile(used_center_value, up_value))]
#   high_name <- used_name[which(used_center_value >= stats::quantile(used_center_value, up_value))]

#   sub_graph <- igraph::induced_subgraph(out$graph, high_name)
#   sub_cpt <- components(sub_graph)
#   cluster <- list()
#   length <- c()
#   pch <- rep(1, length(high))
#   dist_cor <- c()
#   dist_mean <- c()

#   i <- 1
#   while (i <= sub_cpt$no) {
#     this_name <- names(which(sub_cpt$membership == i))
#     this_index <- as.numeric(stringr::str_replace(this_name, "P", ""))
#     pch[which(high %in% this_index)] <- as.character(i)
#     length <- c(length, length(this_index))
#     if (length(this_index) == 1) {
#       this_dist <- dist[used, this_index]
#     } else {
#       this_dist <- apply(dist[used, this_index], 1, mean)
#     }
#     this_cor <- cor(this_dist, center_value[used], method = "spearman")
#     dist_cor <- c(dist_cor, this_cor)

#     dist_mean <- c(dist_mean, mean(this_dist))
#     cluster <- c(cluster, list(this_index))
#     i <- i + 1
#   }

#   tmp <- -dist_cor
#   select <- which(tmp == max(tmp))[1]

#   summit <- cluster[[select]]

#   score <- c()
#   PS <- c()
#   i <- 1
#   while (i <= length(used_name)) {
#     this_name <- used_name[i]
#     this_dist <- dist[
#       which(colnames(dist) == this_name),
#       which(rownames(dist) %in% paste0("P", summit))
#     ]
#     this_dist <- this_dist
#     this_score <- min(this_dist)
#     score <- c(score, this_score)
#     PS <- c(PS, this_score)
#     i <- i + 1
#   }
#   score <- max(score) - score

#   value <- score
#   value_norm <- normalization(value)
#   colors <- vector_vcol(value_norm, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

#   out$colors <- rep("grey70", nrow(out$vector))
#   out$ORIG.COL <- rep("grey70", nrow(out$vector))
#   out$p_score <- rep(0, nrow(out$vector))
#   out$P.PS <- rep(NA, nrow(out$vector))
#   i <- 1
#   while (i <= length(used)) {
#     this_index <- index_list[[used[i]]]
#     out$colors[this_index] <- colors[i]
#     out$ORIG.COL[this_index] <- out$ORIG.CENTER.COL[used][i]
#     out$p_score[this_index] <- score[i]
#     out$P.PS[this_index] <- PS[i]
#     i <- i + 1
#   }

#   if (plot == TRUE) {
#     graphics::plot(
#       out$vector,
#       col = out$ORIG.COL,
#       pch = 16,
#       cex = 0.5
#     )
#     graphics::text(
#       center_vector[high, 1],
#       center_vector[high, 2],
#       labels = pch,
#       cex = 1,
#       pos = 2
#     )
#     graphics::points(
#       center_vector[high, 1],
#       center_vector[high, 2],
#       col = "black",
#       pch = 16,
#       cex = 1
#     )
#     graphics::points(
#       center_vector[summit, 1],
#       center_vector[summit, 2],
#       col = "black",
#       pch = 16,
#       cex = 1.5
#     )
#     graphics::points(
#       center_vector[summit, 1],
#       center_vector[summit, 2],
#       col = "red",
#       pch = 16,
#       cex = 1
#     )
#   }

#   out$score <- score
#   out$summit <- summit
#   out$cluster <- cluster
#   out$length <- length
#   out$pch <- pch
#   out$dist_cor <- dist_cor
#   out$PS <- PS

#   return(out)
# }
