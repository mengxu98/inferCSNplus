vector.SeuratSelect <- function(object) {
  p <- Seurat::DimPlot(
    object,
    reduction = "umap",
    pt.size = 0.5
  )
  used_cells <- Seurat::CellSelector(plot = p)

  return(used_cells)
}

vector.SeuratAddMetaByCell <- function(object, used_cells) {
  SELECT <- rep("NO", ncol(object))
  SELECT[which(colnames(object) %in% used_cells)] <- "YES"
  object@meta.data$select <- SELECT
  return(object)
}

.normX <- function(x) {
  y <- (x - min(x)) / (max(x) - min(x))
  return(y)
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
  #################################
  N.VALUE <- .normX(VALUE)
  COL <- vector.vcol(N.VALUE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))
  graphics::plot(VEC, col = COL, pch = 16, cex = 0.5, xlab = "", ylab = "", axes = FALSE, asp = 1)
  graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NA, border = "black")
  return(OUT)
}

vector.calValue <- function(PCA) {
  OUT <- list()
  PCA.RC <- apply(apply(PCA, 2, rank), 2, .normX)
  PCA.RC <- abs(PCA.RC - 0.5)
  VALUE <- apply(PCA.RC, 1, mean)
  # VALUE=rank(VALUE)
  ###########################
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


#################################################################
# Manually do something

vector.selectPoint <- function(VEC, CEX = 0.5) {
  points <- gatepoints::fhs(VEC, pch = 16, col = "red3", cex = CEX, mark = TRUE)
  return(points)
}

vector.selectCenter <- function(OUT, SHOW = TRUE) {
  DIST <- OUT$DIST
  USED <- OUT$USED
  USED_NAME <- OUT$USED_NAME
  # CENTER_VALUE <- OUT$CENTER_VALUE
  CENTER_VEC <- OUT$CENTER_VEC
  # CENTER_INDEX <- OUT$CENTER_INDEX
  INDEX_LIST <- OUT$INDEX_LIST

  graphics::plot(OUT$VEC, col = "grey80", pch = 16, cex = 0.5, asp = 1)
  graphics::points(OUT$CENTER_VEC[USED, ], col = OUT$ORIG.CENTER.COL[USED], pch = 16, cex = 1)
  SELECT_NAME <- vector.selectPoint(OUT$CENTER_VEC[USED, ], CEX = 1)
  SUMMIT <- USED[as.numeric(SELECT_NAME)]
  # plot(OUT$CENTER_VEC[USED,])
  # graphics::points(OUT$CENTER_VEC[CLUSTER[[9]],],pch=16,col='red')

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
  # plot(OUT$VEC, col='grey70',pch=16)
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

  if (SHOW == TRUE) {
    graphics::plot(OUT$VEC, col = OUT$COL, pch = 16, cex = 0.5, asp = 1)
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
  # A_COL=OUT$A_COL
  VEC <- OUT$VEC

  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16, asp = 1)
  i <- 1
  while (i <= length(A_LENGTH)) {
    graphics::arrows(
      x0 = A1_VEC[i, 1], y0 = A1_VEC[i, 2],
      x1 = A2_VEC[i, 1], y1 = A2_VEC[i, 2],
      lwd = 2, length = A_LENGTH[i],
      col = "black"
      # col=A_COL[i]
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

  # CENTER_LIST <- OUT$CENTER_LIST
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
  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16, asp = 1)
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

  # graphics::points(x=VEC[,1],y=VEC[,2], col='grey70',cex=0.5, pch=16)
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
  # CENTER_LIST <- OUT$CENTER_LIST
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
  graphics::plot(x = VEC[, 1], y = VEC[, 2], col = COL, cex = 0.5, pch = 16, asp = 1)
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
  # graphics::points(x=VEC[,1],y=VEC[,2], col='grey70',cex=0.5, pch=16)
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
  PCA.OUT <- gmodels::fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
  PPCA <- PCA.OUT$x
  return(PPCA)
}

vector.SeuratPCA <- function(object, CUT = 0.5) {
  D <- as.matrix(object@assays$RNA@scale.data)
  PCA.OUT <- gmodels::fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
  EXP <- (cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2))
  N <- min(which(cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2) > CUT))
  PCA.OUT$EXP <- EXP
  PCA.OUT$CUT <- CUT
  PCA.OUT$N <- N
  # N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > 0.7))

  return(PCA.OUT)
}

vector.SeuratRandomPCA <- function(object, RN = 1000, CUT = 0.5) {
  D <- as.matrix(object@assays$RNA@scale.data)
  R_INDEX <- sample(1:ncol(D), RN, replace = TRUE)

  RD <- D[, R_INDEX]
  PCA.OUT <- gmodels::fast.prcomp(t(RD), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
  EXP <- (cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2))
  N <- min(which(cumsum(PCA.OUT$sdev^2) / sum(PCA.OUT$sdev^2) > CUT))

  PRED.PCA <- t(D) %*% PCA.OUT$rotation

  PCA.OUT$EXP <- EXP
  PCA.OUT$CUT <- CUT
  PCA.OUT$N <- N
  PCA.OUT$RN <- RN
  PCA.OUT$R_INDEX <- R_INDEX
  PCA.OUT$PRED.PCA <- PRED.PCA
  # N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > 0.7))

  return(PCA.OUT)
}

vector.medCurv <- function(PCA, MAX = 1000) {
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

vector.getValueNew <- function(OUT, PCA, SHOW = TRUE) {
  VALUE.OUT <- vector.calValueNew(PCA)
  OUT$VALUE <- VALUE.OUT$VALUE
  OUT$PCA <- PCA
  OUT$PCA.RC <- VALUE.OUT$PCA.RC
  if (SHOW == TRUE) {
    vector.showValue(OUT)
  }

  return(OUT)
}


# 20201030
vector.autoCenterCor <- function(OUT, UP = 0.9, SHOW = TRUE) {
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

  # SELECT=which(DIST_COR==min(DIST_COR))[1]
  # SELECT=which(rank(DIST_MEAN)*rank(DIST_COR) == min(rank(DIST_MEAN)*rank(DIST_COR)) )

  # TMP=10^LENGTH - DIST_COR
  TMP <- -DIST_COR
  # TMP=rank(-DIST_COR) * rank(LENGTH )
  SELECT <- which(TMP == max(TMP))[1]
  # SELECT=which(LENGTH==max(LENGTH))[1]
  # SELECT=which( rank(-DIST_COR) * rank(LENGTH) == max( rank(-DIST_COR) * rank(LENGTH)  ) )[1]

  SUMMIT <- CLUSTER[[SELECT]]
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
    this_dist <- this_dist
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

  if (SHOW == TRUE) {
    graphics::plot(OUT$VEC, col = OUT$ORIG.COL, pch = 16, cex = 0.5, asp = 1)
    graphics::text(CENTER_VEC[HIGH, 1], CENTER_VEC[HIGH, 2], labels = PCH, cex = 1, pos = 2)
    graphics::points(CENTER_VEC[HIGH, 1], CENTER_VEC[HIGH, 2], col = "black", pch = 16, cex = 1)
    graphics::points(CENTER_VEC[SUMMIT, 1], CENTER_VEC[SUMMIT, 2], col = "black", pch = 16, cex = 1.5)
    graphics::points(CENTER_VEC[SUMMIT, 1], CENTER_VEC[SUMMIT, 2], col = "red", pch = 16, cex = 1)
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
