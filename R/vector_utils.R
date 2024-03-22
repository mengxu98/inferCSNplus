#' rankPCA
#'
#' @param PCA input
#'
#' @return output
#' @export
vector.rankPCA <- function(PCA) {
  R.PCA <- apply(PCA, 2, rank)
  PCA.OUT <- gmodels::fast.prcomp(
    R.PCA,
    retx = TRUE,
    center = TRUE,
    scale. = TRUE,
    tol = NULL
  )
  N.PCA <- PCA.OUT$x

  return(N.PCA)
}

#' Title
#'
#' @param VEC VEC
#' @param N N
#' @param plot plot
#' @param COL COL
#'
#' @return output
#' @export
vector.buildGrid <- function(
    VEC,
    N = 30,
    plot = TRUE,
    COL = "grey70") {
  VEC.E <- VEC
  delta <- 0.000001
  N <- N

  if (plot == TRUE) {
    graphics::plot(
      VEC.E,
      col = COL,
      pch = 16,
      cex = 0.2,
      xlab = "",
      ylab = "",
      axes = FALSE)
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black")
  }

  X <- VEC[, 1]
  Y <- VEC[, 2]

  this_step_x <- (max(X) - min(X) - delta) / N
  this_step_y <- (max(Y) - min(Y) - delta) / N
  NUM_CUT <- 1

  INDEX_LIST <- list()
  CENTER_LIST <- list()

  this_x <- min(X)
  while (this_x < max(X)) {
    this_y <- min(Y)
    while (this_y < max(Y)) {
      this_in_index <- which(VEC.E[, 1] >= this_x & VEC.E[, 1] < this_x + this_step_x &
                               VEC.E[, 2] >= this_y & VEC.E[, 2] < this_y + this_step_y)
      this_center <- c(this_x + this_step_x / 2, this_y + this_step_y / 2)
      if (length(this_in_index) >= NUM_CUT) {
        INDEX_LIST <- c(INDEX_LIST, list(this_in_index))
        CENTER_LIST <- c(CENTER_LIST, list(this_center))
        if (plot == TRUE) {
          graphics::points(this_center[1], this_center[2], col = "black", pch = 16, cex = 0.5)
        }
      }
      this_y <- this_y + this_step_y
    }
    this_x <- this_x + this_step_x
  }

  OUT <- list()
  OUT$VEC <- VEC.E
  OUT$INDEX_LIST <- INDEX_LIST
  OUT$CENTER_LIST <- CENTER_LIST
  OUT$this_step_x <- this_step_x
  OUT$this_step_y <- this_step_y

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param CUT CUT
#' @param plot plot
#' @param COL COL
#'
#' @return output
#' @export
vector.buildNet <- function(
    OUT,
    CUT = 1,
    plot = TRUE,
    COL = "grey70") {
  VEC <- OUT$VEC
  INDEX_LIST <- OUT$INDEX_LIST
  CENTER_LIST <- OUT$CENTER_LIST
  this_step_x <- OUT$this_step_x
  this_step_y <- OUT$this_step_y
  delta <- 0.00001

  if (plot == TRUE) {
    graphics::plot(
      VEC,
      col = COL,
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

  CENTER_VEC <- c()
  i <- 1
  while (i <= length(CENTER_LIST)) {
    CENTER_VEC <- cbind(CENTER_VEC, CENTER_LIST[[i]])
    i <- i + 1
  }
  CENTER_VEC <- t(CENTER_VEC)
  OUT$CENTER_VEC <- CENTER_VEC

  p1 <- c()
  p2 <- c()

  CNUM <- length(CENTER_LIST)
  i <- 1
  while (i <= CNUM) {
    this_p1_loc <- CENTER_LIST[[i]]
    this_p1 <- paste0("P", as.character(i))

    used_j <- which(
      (abs(CENTER_VEC[, 1] - this_p1_loc[1]) <= this_step_x + delta) &
        (abs(CENTER_VEC[, 2] - this_p1_loc[2]) <= this_step_y + delta)
    )

    for (j in used_j) {
      this_p2 <- paste0("P", as.character(j))

      this_p2_loc <- CENTER_LIST[[j]]

      if (length(INDEX_LIST[[i]]) >= CUT &&
          length(INDEX_LIST[[j]]) >= CUT &&
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

  OUT$p1 <- p1
  OUT$p2 <- p2
  NET <- cbind(p1, p2)
  g <- igraph::make_graph(t(NET), directed = FALSE)
  ALLNAME <- paste0("P", 1:CNUM)
  ADD <- ALLNAME[which(!ALLNAME %in% igraph::as_ids(igraph::V(g)))]
  g <- g + ADD

  DIST <- igraph::distances(
    g,
    v = igraph::V(g),
    to = igraph::V(g),
    mode = c("all")
  )

  DIST.NUM <- as.numeric(
    stringr::str_replace(colnames(DIST), "P", "")
  )
  DIST <- DIST[order(DIST.NUM), order(DIST.NUM)]

  CPT <- igraph::components(g)
  MAXC <- which(CPT$csize == max(CPT$csize))[1]

  USED_NAME <- names(which(CPT$membership == MAXC))
  USED <- as.numeric(
    stringr::str_replace(USED_NAME, "P", "")
  )

  USED_NAME <- USED_NAME[order(USED)]
  USED <- USED[order(USED)]
  if (plot == TRUE) {
    graphics::points(OUT$CENTER_VEC[USED, ], col = "red", pch = 16, cex = 0.5)
  }
  USED_INDEX <- c()
  i <- 1
  while (i <= length(USED)) {
    USED_INDEX <- c(USED_INDEX, INDEX_LIST[[USED[i]]])
    i <- i + 1
  }

  OUT$GRAPH <- g
  OUT$DIST <- DIST
  OUT$USED <- USED
  OUT$USED_NAME <- USED_NAME
  OUT$USED_INDEX <- USED_INDEX

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param PCA PCA
#' @param plot plot
#'
#' @return output
#' @export
vector.getValue <- function(OUT, PCA, plot = TRUE) {
  VALUE.OUT <- vector.calValue(PCA)

  OUT$VALUE <- VALUE.OUT$VALUE
  OUT$PCA <- PCA
  OUT$PCA.RC <- VALUE.OUT$PCA.RC

  if (plot == TRUE) {
    vector.showValue(OUT)
  }

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param plot plot
#'
#' @return output
#' @export
vector.gridValue <- function(OUT, plot = TRUE) {
  INDEX_LIST <- OUT$INDEX_LIST
  VALUE <- OUT$VALUE
  plot <- plot
  USED <- OUT$USED

  CENTER_VALUE <- c()
  i <- 1
  while (i <= length(INDEX_LIST)) {
    this_value <- mean(VALUE[INDEX_LIST[[i]]])
    CENTER_VALUE <- c(CENTER_VALUE, this_value)
    i <- i + 1
  }
  # CENTER_VEC <- OUT$CENTER_VEC

  N.VALUE <- (VALUE - min(VALUE)) / (max(VALUE) - min(VALUE))
  VALUE.COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  VALUE <- CENTER_VALUE
  N.VALUE <- (VALUE - min(VALUE)) / (max(VALUE) - min(VALUE))
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  if (plot == TRUE) {
    graphics::plot(
      OUT$VEC,
      col = "grey80",
      pch = 16,
      cex = 0.5,
      xlab = "",
      ylab = "",
      axes = FALSE)
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black")
    # graphics::points(OUT$CENTER_VEC,col=COL,pch=16)
    graphics::points(OUT$CENTER_VEC[USED, ], col = COL[USED], pch = 15, cex = 1.5)
  }
  OUT$CENTER_VALUE <- CENTER_VALUE
  OUT$ORIG.CENTER.COL <- COL
  OUT$VALUE.COL <- VALUE.COL

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param UP UP
#' @param plot plot
#'
#' @return output
#' @export
vector.autoCenter <- function(OUT, UP = 0.9, plot = TRUE) {
  DIST <- OUT$DIST
  USED <- OUT$USED
  USED_NAME <- OUT$USED_NAME
  CENTER_VALUE <- OUT$CENTER_VALUE
  CENTER_VEC <- OUT$CENTER_VEC
  # CENTER_INDEX <- OUT$CENTER_INDEX
  INDEX_LIST <- OUT$INDEX_LIST

  USED_CENTER_VALUE <- CENTER_VALUE[USED]
  # USED_CENTER_VEC <- CENTER_VEC[USED, ]

  HIGH <- USED[which(USED_CENTER_VALUE >= stats::quantile(USED_CENTER_VALUE, UP))]
  HIGH_NAME <- USED_NAME[which(USED_CENTER_VALUE >= stats::quantile(USED_CENTER_VALUE, UP))]
  # LOW_NAME <- USED_NAME[which(USED_CENTER_VALUE < stats::quantile(USED_CENTER_VALUE, UP))]
  # plot(OUT$CENTER_VEC[USED,])
  # graphics::points(CENTER_VEC[HIGH,], col='red',pch=16)
  # plot(CENTER_VEC[HIGH,], col='red')

  SUB <- igraph::induced_subgraph(OUT$GRAPH, HIGH_NAME)
  SUB_CPT <- igraph::components(SUB)
  CLUSTER <- list()
  LENGTH <- c()
  PCH <- rep(1, length(HIGH))
  DIST_COR <- c()
  DIST_MEAN <- c()

  i <- 1
  while (i <= SUB_CPT$no) {
    this_name <- names(which(SUB_CPT$membership == i))
    this_index <- as.numeric(stringr::str_replace(this_name, "P", ""))
    PCH[which(HIGH %in% this_index)] <- as.character(i)
    LENGTH <- c(LENGTH, length(this_index))
    if (length(this_index) == 1) {
      this_dist <- DIST[USED, this_index]
    } else {
      this_dist <- apply(DIST[USED, this_index], 1, mean)
    }
    this_cor <- cor(this_dist, CENTER_VALUE[USED], method = "spearman")
    DIST_COR <- c(DIST_COR, this_cor)

    DIST_MEAN <- c(DIST_MEAN, mean(this_dist))
    CLUSTER <- c(CLUSTER, list(this_index))
    i <- i + 1
  }

  # SELECT=which(DIST_COR==min(DIST_COR))[1]
  # SELECT=which(rank(DIST_MEAN)*rank(DIST_COR) == min(rank(DIST_MEAN)*rank(DIST_COR)) )

  TMP <- 10^LENGTH - DIST_COR
  # TMP=rank(-DIST_COR) * rank(LENGTH )
  SELECT <- which(TMP == max(TMP))[1]
  # SELECT=which(LENGTH==max(LENGTH))[1]
  # SELECT=which( rank(-DIST_COR) * rank(LENGTH) == max( rank(-DIST_COR) * rank(LENGTH)  ) )[1]

  SUMMIT <- CLUSTER[[SELECT]]
  # print(SELECT)
  # plot(OUT$CENTER_VEC[USED,])
  # graphics::points(OUT$CENTER_VEC[CLUSTER[[9]],],pch=16,col='red')

  SCORE <- c()
  PS <- c()
  i <- 1
  while (i <= length(USED_NAME)) {
    this_name <- USED_NAME[i]
    this_dist <- DIST[
      which(colnames(DIST) == this_name),
      which(rownames(DIST) %in% paste0("P", SUMMIT))
    ]
    # this_value=CENTER_VALUE[USED]
    this_score <- min(this_dist) # sum(rank(-this_dist) * rank(this_value))
    # this_cor=cor(-this_value, this_dist)#,method='spearman')
    SCORE <- c(SCORE, this_score)
    PS <- c(PS, this_score)
    i <- i + 1
  }
  SCORE <- max(SCORE) - SCORE

  VALUE <- SCORE
  # plot(OUT$VEC, col='grey70',pch=16)
  N.VALUE <- (VALUE - min(VALUE)) / (max(VALUE) - min(VALUE))
  # COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) # too strong
  # COL=vector.vcol(N.VALUE,c(0,0.5,1),c('#009FFF','#FFF200','#ffdde1')) # too weak
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ee9ca7"))

  OUT$COL <- rep("grey70", nrow(OUT$VEC))
  OUT$ORIG.COL <- rep("grey70", nrow(OUT$VEC))
  OUT$P.SCORE <- rep(0, nrow(OUT$VEC))
  OUT$P.PS <- rep(NA, nrow(OUT$VEC))
  i <- 1
  while (i <= length(USED)) {
    this_index <- INDEX_LIST[[USED[i]]]
    OUT$COL[this_index] <- COL[i]
    OUT$ORIG.COL[this_index] <- OUT$ORIG.CENTER.COL[USED][i]
    OUT$P.SCORE[this_index] <- SCORE[i]
    OUT$P.PS[this_index] <- PS[i]
    i <- i + 1
  }

  if (plot == TRUE) {
    graphics::plot(
      OUT$VEC,
      col = OUT$ORIG.COL,
      pch = 16,
      cex = 0.5,
      xlab = "",
      ylab = "",
      axes = FALSE)
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black")
    graphics::text(
      CENTER_VEC[HIGH, 1],
      CENTER_VEC[HIGH, 2],
      labels = PCH,
      cex = 1,
      pos = 2)
    graphics::points(
      CENTER_VEC[HIGH, 1],
      CENTER_VEC[HIGH, 2],
      col = "black",
      pch = 16,
      cex = 1)
    graphics::points(
      CENTER_VEC[SUMMIT, 1],
      CENTER_VEC[SUMMIT, 2],
      col = "black",
      pch = 16,
      cex = 1.5)
    graphics::points(
      CENTER_VEC[SUMMIT, 1],
      CENTER_VEC[SUMMIT, 2],
      col = "red",
      pch = 16,
      cex = 1)
  }

  OUT$SCORE <- SCORE
  OUT$SUMMIT <- SUMMIT
  OUT$CLUSTER <- CLUSTER
  OUT$LENGTH <- LENGTH
  OUT$PCH <- PCH
  OUT$DIST_COR <- DIST_COR
  OUT$PS <- PS
  # OUT$DIST_MEAN=DIST_MEAN

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param P P
#' @param plot plot
#' @param COL COL
#' @param OL OL
#' @param arrow_length2 arrow_length2
#' @param CEX CEX
#' @param arrow_width arrow_width
#' @param BD BD
#' @param arrow_color arrow_color
#' @param plot.SUMMIT plot.SUMMIT
#'
#' @return output
#' @export
vector.drawArrow <- function(
    OUT,
    P = 0.9,
    plot = TRUE,
    COL = "grey70",
    OL = 1.5,
    arrow_length2 = 70,
    CEX = 0.5,
    arrow_width = 1,
    BD = TRUE,
    arrow_color = "grey30",
    plot.SUMMIT = TRUE) {
  USED <- OUT$USED
  DIST <- OUT$DIST
  ALL_VEC <- OUT$VEC
  USED_NAME <- OUT$USED_NAME
  USED_CENTER_VEC <- OUT$CENTER_VEC[USED, ]
  USED_DIST <- OUT$DIST[which(rownames(DIST) %in% USED_NAME), which(rownames(DIST) %in% USED_NAME)]
  SCORE <- OUT$SCORE

  one <- min(stats::dist(USED_CENTER_VEC)) * OL

  DIV <- 1 / P

  if (plot == TRUE) {
    if (BD == TRUE) {
      graphics::plot(
        ALL_VEC,
        col = COL,
        pch = 16,
        cex = CEX,
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
        col = COL,
        pch = 16,
        cex = CEX,
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
  N.SCORE <- normalization(SCORE)
  SCORE.COL <- vector.vcol(
    N.SCORE,
    c(0, 0.5, 1),
    c("#009FFF", "#FFF200", "#ec2F4B")
  )

  A1_VEC <- c()
  A2_VEC <- c()
  A_LENGTH <- c()

  i <- 1
  while (i <= length(USED)) {
    this_p1_loc <- USED_CENTER_VEC[i, ]

    vector_list <- cbind(USED_CENTER_VEC[, 1] - this_p1_loc[1], USED_CENTER_VEC[, 2] - this_p1_loc[2])
    vector_list_norm <- t(apply(vector_list, 1, .norm_one, one))

    vector_weight_1 <- DIV^-(rank(USED_DIST[i, ]) - 1) # equals to P ^ (rank(USED_DIST[i,])-1)
    vector_weight_2 <- SCORE[i] - SCORE # PS - PS[i]; "PS_i - PS_j" in  the paper; In paper, PS[i] is "PS_j".

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
    X1 <- min(OUT$CENTER_VEC[OUT$SUMMIT, 1]) - one / 10
    X2 <- max(OUT$CENTER_VEC[OUT$SUMMIT, 1]) + one / 10
    Y1 <- min(OUT$CENTER_VEC[OUT$SUMMIT, 2]) - one / 10
    Y2 <- max(OUT$CENTER_VEC[OUT$SUMMIT, 2]) + one / 10

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
  OUT$A1_VEC <- A1_VEC
  OUT$A2_VEC <- A2_VEC
  OUT$A_LENGTH <- A_LENGTH
  OUT$A_COL <- SCORE.COL

  return(OUT)
}

# Other functions
vector.select.Seurat <- function(object) {
  p <- Seurat::DimPlot(
    object,
    reduction = "umap",
    pt.size = 0.5
  )
  used_cells <- Seurat::CellSelector(plot = p)

  return(used_cells)
}

vector.AddMetaByCell.Seurat <- function(
    object,
    used_cells) {
  SELECT <- rep("NO", ncol(object))
  SELECT[which(colnames(object) %in% used_cells)] <- "YES"
  object@meta.data$select <- SELECT
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

  PCA.OUT <- switch(
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

  EXP <- (cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2))
  N <- min(which(cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2) > CUT))

  if (random) {
    PRED.PCA <- t(D_raw) %*% PCA.OUT$rotation
  }

  PCA.OUT$EXP <- EXP
  PCA.OUT$CUT <- CUT
  PCA.OUT$N <- N
  if (random) {
    PCA.OUT$RN <- RN
    PCA.OUT$R_INDEX <- R_INDEX
    PCA.OUT$PRED.PCA <- PRED.PCA
  }

  return(PCA.OUT)
}

vector.lcol <- function(tag) {
  tag <- as.factor(tag)
  my_color_palette <- scales::hue_pal()(length(unique(tag)))
  COL <- my_color_palette[tag]
  return(COL)
}

vector.vcol <- function(VALUE, CV, CN) {
  CRF <- circlize::colorRamp2(CV, CN)
  COL <- CRF(VALUE)
  return(COL)
}

vector.showValue <- function(OUT) {
  VEC <- OUT$VEC
  VALUE <- OUT$VALUE
  N.VALUE <- normalization(VALUE)
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))
  graphics::plot(VEC, col = COL, pch = 16, cex = 0.5, xlab = "", ylab = "", axes = FALSE)
  graphics::rect(
    graphics::par("usr")[1],
    graphics::par("usr")[3],
    graphics::par("usr")[2],
    graphics::par("usr")[4],
    col = NA,
    border = "black")
  return(OUT)
}

vector.calValue <- function(PCA) {
  OUT <- list()
  PCA.RC <- apply(apply(PCA, 2, rank), 2, normalization)
  PCA.RC <- abs(PCA.RC - 0.5)
  VALUE <- apply(PCA.RC, 1, mean)
  OUT$VALUE <- VALUE
  OUT$PCA.RC <- PCA.RC
  return(OUT)
}

vector.nonCenter <- function(OUT) {
  OUT$SCORE <- OUT$CENTER_VALUE[OUT$USED]
  OUT$COL <- OUT$VALUE.COL
  return(OUT)
}

.norm_one <- function(x, one = 1) {
  if (stats::var(x) != 0) {
    x <- x / sqrt(sum(x^2)) * one
  }
  return(x)
}

# Manually do something
vector.selectPoint <- function(VEC, CEX = 0.5) {
  points <- gatepoints::fhs(VEC, pch = 16, col = "red3", cex = CEX, mark = TRUE)
  return(points)
}

vector.selectCenter <- function(OUT, plot = TRUE) {
  DIST <- OUT$DIST
  USED <- OUT$USED
  USED_NAME <- OUT$USED_NAME
  CENTER_VEC <- OUT$CENTER_VEC
  INDEX_LIST <- OUT$INDEX_LIST

  graphics::plot(OUT$VEC, col = "grey80", pch = 16, cex = 0.5)
  graphics::points(OUT$CENTER_VEC[USED, ], col = OUT$ORIG.CENTER.COL[USED], pch = 16, cex = 1)
  SELECT_NAME <- vector.selectPoint(OUT$CENTER_VEC[USED, ], CEX = 1)
  SUMMIT <- USED[as.numeric(SELECT_NAME)]

  SCORE <- c()
  i <- 1
  while (i <= length(USED_NAME)) {
    this_name <- USED_NAME[i]
    this_dist <- DIST[which(colnames(DIST) == this_name), which(rownames(DIST) %in% paste0("P", SUMMIT))]
    this_dist <- this_dist
    # this_value=CENTER_VALUE[USED]
    this_score <- min(this_dist) # sum(rank(-this_dist) * rank(this_value))
    # this_cor=cor(-this_value, this_dist)#,method='spearman')
    SCORE <- c(SCORE, this_score)
    i <- i + 1
  }

  PS <- SCORE

  SCORE <- max(SCORE) - SCORE

  VALUE <- SCORE
  N.VALUE <- (VALUE - min(VALUE)) / (max(VALUE) - min(VALUE))
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  OUT$COL <- rep("grey70", nrow(OUT$VEC))
  OUT$ORIG.COL <- rep("grey70", nrow(OUT$VEC))
  OUT$P.SCORE <- rep(0, nrow(OUT$VEC))
  OUT$P.PS <- rep(NA, nrow(OUT$VEC))
  i <- 1
  while (i <= length(USED)) {
    this_index <- INDEX_LIST[[USED[i]]]
    OUT$COL[this_index] <- COL[i]
    OUT$ORIG.COL[this_index] <- OUT$ORIG.CENTER.COL[USED][i]
    OUT$P.SCORE[this_index] <- SCORE[i]
    OUT$P.PS[this_index] <- PS[i]
    i <- i + 1
  }

  if (plot == TRUE) {
    graphics::plot(OUT$VEC, col = OUT$COL, pch = 16, cex = 0.5)
    graphics::points(CENTER_VEC[SUMMIT, ], col = "black", pch = 16, cex = 1.5)
    graphics::points(CENTER_VEC[SUMMIT, ], col = "red", pch = 16, cex = 1)
  }

  OUT$SCORE <- SCORE
  OUT$PS <- PS
  OUT$SUMMIT <- SUMMIT

  return(OUT)
}

vector.reDrawArrow <- function(OUT, COL = "grey70") {
  A1_VEC <- OUT$A1_VEC
  A2_VEC <- OUT$A2_VEC
  A_LENGTH <- OUT$A_LENGTH
  VEC <- OUT$VEC

  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16)
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

  return(OUT)
}


vector.selectRegion <- function(OUT) {
  P.SCORE <- OUT$P.SCORE
  A1_VEC <- OUT$A1_VEC
  A2_VEC <- OUT$A2_VEC
  A_LENGTH <- OUT$A_LENGTH
  VEC <- OUT$VEC
  P.PS <- OUT$P.PS

  INDEX_LIST <- OUT$INDEX_LIST
  USED <- OUT$USED
  SELECT_NAME <- vector.selectPoint(VEC, CEX = 0.1)
  SELECT_INDEX <- which(rownames(VEC) %in% SELECT_NAME)
  SELECT_INDEX <- SELECT_INDEX[which(SELECT_INDEX %in% OUT$USED_INDEX)]

  A_USED <- c()
  i <- 1
  while (i <= length(USED)) {
    if (length(which(INDEX_LIST[[USED[i]]] %in% SELECT_INDEX)) > 0) {
      A_USED <- c(A_USED, i)
    }

    i <- i + 1
  }

  # Draw new
  COL <- OUT$COL
  COL[which(!rownames(VEC) %in% SELECT_NAME)] <- "grey70"
  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16)
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

  OUT$A_USED <- A_USED
  OUT$SELECT_NAME <- SELECT_NAME
  OUT$SELECT_INDEX <- SELECT_INDEX
  OUT$SELECT_SCORE <- P.SCORE[SELECT_INDEX]
  OUT$SELECT_PS <- P.PS[SELECT_INDEX]

  return(OUT)
}

vector.reDrawRegion <- function(OUT) {
  P.SCORE <- OUT$P.SCORE
  A1_VEC <- OUT$A1_VEC
  A2_VEC <- OUT$A2_VEC
  A_LENGTH <- OUT$A_LENGTH
  VEC <- OUT$VEC
  P.PS <- OUT$P.PS
  INDEX_LIST <- OUT$INDEX_LIST
  USED <- OUT$USED
  SELECT_NAME <- OUT$SELECT_NAME # vector.selectPoint(VEC,CEX=0.1)
  SELECT_INDEX <- which(rownames(VEC) %in% SELECT_NAME)

  A_USED <- c()
  i <- 1
  while (i <= length(USED)) {
    if (length(which(INDEX_LIST[[USED[i]]] %in% SELECT_INDEX)) > 0) {
      A_USED <- c(A_USED, i)
    }

    i <- i + 1
  }

  # Draw new
  COL <- OUT$COL
  COL[which(!rownames(VEC) %in% SELECT_NAME)] <- "grey70"
  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16)
  # graphics::points(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)

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
  OUT$A_USED <- A_USED
  OUT$SELECT_NAME <- SELECT_NAME
  OUT$SELECT_INDEX <- SELECT_INDEX
  OUT$SELECT_SCORE <- P.SCORE[SELECT_INDEX]
  OUT$SELECT_PS <- P.PS[SELECT_INDEX]

  return(OUT)
}

vector.RPPCA <- function(PCA) {
  R.PCA <- apply(PCA, 2, rank)
  D <- t(R.PCA)
  PCA.OUT <- gmodels::fast.prcomp(
    t(D),
    retx = TRUE,
    center = FALSE,
    scale. = FALSE,
    tol = NULL)
  PPCA <- PCA.OUT$x
  return(PPCA)
}


vector.medCurv <- function(
    PCA,
    MAX = 1000) {
  MS <- vector.calValue(PCA)$PCA.RC
  COR <- cor(MS)
  MED <- c()
  ALL <- c()
  i <- 1
  while (i <= ncol(COR) & i <= MAX) {
    this_col <- COR[1:i, i]
    ALL <- c(ALL, this_col)
    ALL[which(ALL == 1)] <- NA
    this_med <- median(ALL, na.rm = TRUE)
    MED <- c(MED, this_med)
    if (i %% 100 == 1) {
      print(i)
    }
    i <- i + 1
  }
  MED[1] <- MED[2]
  return(MED)
}

# 2020.01.03
vector.smoothOut <- function(X, Z) {
  Y <- X
  Y[order(Z)] <- X[order(Z)] - stats::smooth(X[order(Z)])
  return(Y)
}

vector.regressOut <- function(X, Z) {
  FIT <- stats::lm(X ~ Z)
  Y <- X - predict(FIT)
  return(Y)
}

vector.removeOut <- function(X) {
  Q1 <- stats::quantile(X, 0.25)
  Q3 <- stats::quantile(X, 0.75)
  IQR <- Q3 - Q1
  LW <- Q1 - 1.5 * IQR
  UP <- Q3 + 1.5 * IQR
  X[which(X > UP)] <- UP
  X[which(X < LW)] <- LW

  return(X)
}

# 2020.1.12
.normNew <- function(x) {
  pos_index <- which(x > 0)
  neg_index <- which(x < 0)
  x[pos_index] <- rank(x[pos_index])
  x[neg_index] <- rank(-x[neg_index])

  return(x)
}

vector.calValueNew <- function(PCA) {
  OUT <- list()
  PCA.RC <- apply(PCA, 2, .normNew)
  VALUE <- apply(PCA.RC, 1, mean)
  OUT$VALUE <- VALUE
  OUT$PCA.RC <- PCA.RC

  return(OUT)
}

vector.getValueNew <- function(OUT, PCA, plot = TRUE) {
  VALUE.OUT <- vector.calValueNew(PCA)
  OUT$VALUE <- VALUE.OUT$VALUE
  OUT$PCA <- PCA
  OUT$PCA.RC <- VALUE.OUT$PCA.RC
  if (plot == TRUE) {
    vector.showValue(OUT)
  }

  return(OUT)
}


# 20201030
vector.autoCenterCor <- function(
    OUT,
    UP = 0.9,
    plot = TRUE) {
  DIST <- OUT$DIST
  USED <- OUT$USED
  USED_NAME <- OUT$USED_NAME
  CENTER_VALUE <- OUT$CENTER_VALUE
  CENTER_VEC <- OUT$CENTER_VEC
  INDEX_LIST <- OUT$INDEX_LIST

  USED_CENTER_VALUE <- CENTER_VALUE[USED]

  HIGH <- USED[which(USED_CENTER_VALUE >= stats::quantile(USED_CENTER_VALUE, UP))]
  HIGH_NAME <- USED_NAME[which(USED_CENTER_VALUE >= stats::quantile(USED_CENTER_VALUE, UP))]

  SUB <- igraph::induced_subgraph(OUT$GRAPH, HIGH_NAME)
  SUB_CPT <- components(SUB)
  CLUSTER <- list()
  LENGTH <- c()
  PCH <- rep(1, length(HIGH))
  DIST_COR <- c()
  DIST_MEAN <- c()

  i <- 1
  while (i <= SUB_CPT$no) {
    this_name <- names(which(SUB_CPT$membership == i))
    this_index <- as.numeric(stringr::str_replace(this_name, "P", ""))
    PCH[which(HIGH %in% this_index)] <- as.character(i)
    LENGTH <- c(LENGTH, length(this_index))
    if (length(this_index) == 1) {
      this_dist <- DIST[USED, this_index]
    } else {
      this_dist <- apply(DIST[USED, this_index], 1, mean)
    }
    this_cor <- cor(this_dist, CENTER_VALUE[USED], method = "spearman")
    DIST_COR <- c(DIST_COR, this_cor)

    DIST_MEAN <- c(DIST_MEAN, mean(this_dist))
    CLUSTER <- c(CLUSTER, list(this_index))
    i <- i + 1
  }

  TMP <- -DIST_COR
  SELECT <- which(TMP == max(TMP))[1]

  SUMMIT <- CLUSTER[[SELECT]]

  SCORE <- c()
  PS <- c()
  i <- 1
  while (i <= length(USED_NAME)) {
    this_name <- USED_NAME[i]
    this_dist <- DIST[
      which(colnames(DIST) == this_name),
      which(rownames(DIST) %in% paste0("P", SUMMIT))
    ]
    this_dist <- this_dist
    this_score <- min(this_dist)
    SCORE <- c(SCORE, this_score)
    PS <- c(PS, this_score)
    i <- i + 1
  }
  SCORE <- max(SCORE) - SCORE

  VALUE <- SCORE
  N.VALUE <- (VALUE - min(VALUE)) / (max(VALUE) - min(VALUE))
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  OUT$COL <- rep("grey70", nrow(OUT$VEC))
  OUT$ORIG.COL <- rep("grey70", nrow(OUT$VEC))
  OUT$P.SCORE <- rep(0, nrow(OUT$VEC))
  OUT$P.PS <- rep(NA, nrow(OUT$VEC))
  i <- 1
  while (i <= length(USED)) {
    this_index <- INDEX_LIST[[USED[i]]]
    OUT$COL[this_index] <- COL[i]
    OUT$ORIG.COL[this_index] <- OUT$ORIG.CENTER.COL[USED][i]
    OUT$P.SCORE[this_index] <- SCORE[i]
    OUT$P.PS[this_index] <- PS[i]
    i <- i + 1
  }

  if (plot == TRUE) {
    graphics::plot(
      OUT$VEC,
      col = OUT$ORIG.COL,
      pch = 16,
      cex = 0.5)
    graphics::text(
      CENTER_VEC[HIGH, 1],
      CENTER_VEC[HIGH, 2],
      labels = PCH,
      cex = 1,
      pos = 2)
    graphics::points(
      CENTER_VEC[HIGH, 1],
      CENTER_VEC[HIGH, 2],
      col = "black",
      pch = 16,
      cex = 1)
    graphics::points(
      CENTER_VEC[SUMMIT, 1],
      CENTER_VEC[SUMMIT, 2],
      col = "black",
      pch = 16,
      cex = 1.5)
    graphics::points(
      CENTER_VEC[SUMMIT, 1],
      CENTER_VEC[SUMMIT, 2],
      col = "red",
      pch = 16,
      cex = 1)
  }

  OUT$SCORE <- SCORE
  OUT$SUMMIT <- SUMMIT
  OUT$CLUSTER <- CLUSTER
  OUT$LENGTH <- LENGTH
  OUT$PCH <- PCH
  OUT$DIST_COR <- DIST_COR
  OUT$PS <- PS

  return(OUT)
}
