#' vector.main
#'
#' @description
#' Raw source code form: \url{https://github.com/jumphone/Vector}.
#' Author: Feng Zhang.
#' Paper: Unsupervised Inference of Developmental Directions for Single Cells Using VECTOR,
#' please reference: \url{https://doi.org/10.1016/j.celrep.2020.108069}.
#'
#' @param embedding_information The information of PCA and UMAP.
#' @param plot Logical value.
#'
#' @return vector_result
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' embedding_information <- get.embedding(example_matrix)
#' vector_result <- vector.main(embedding_information)
vector.main <- function(
    embedding_information,
    plot = TRUE) {
  PCA <- embedding_information[["PCA"]]
  VEC <- embedding_information[["UMAP"]]

  graphics::par(mfrow = c(2, 3), mai = c(0.2, 0.2, 0.1, 0.1))

  PCA <- vector.rankPCA(PCA)

  # Define pixel
  OUT <- vector.buildGrid(VEC, N = 30, SHOW = plot)

  # Build network
  OUT <- vector.buildNet(OUT, CUT = 1, SHOW = plot)

  # Calculate Quantile Polarization (QP) score
  OUT <- vector.getValue(OUT, PCA, SHOW = plot)

  # Get pixel's QP score
  OUT <- vector.gridValue(OUT, SHOW = plot)

  # Find starting point
  OUT <- vector.autoCenter(OUT, UP = 0.9, SHOW = plot)

  # Infer vector
  OUT <- vector.drawArrow(
    OUT,
    P = 0.9,
    SHOW = plot,
    COL = OUT$COL,
    SHOW.SUMMIT = plot
  )

  # OUT$P.PS : Peseudotime Score (PS) of each cell

  graphics::par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

  graphics::mtext("UMAP1", side = 1)
  graphics::mtext("UMAP2", side = 2)

  return(OUT)
}

#' rankPCA
#'
#' @param PCA input
#'
#' @return output
#' @export
vector.rankPCA <- function(PCA) {
  PCA <- PCA # ncol=num(PC)
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
#' @param SHOW SHOW
#' @param COL COL
#'
#' @return output
#' @export
vector.buildGrid <- function(
    VEC,
    N = 30,
    SHOW = TRUE,
    COL = "grey70") {
  VEC.E <- VEC
  delta <- 0.000001
  N <- N

  if (SHOW == TRUE) {
    graphics::plot(VEC.E, col = COL, pch = 16, cex = 0.2, xlab = "", ylab = "", axes = FALSE, asp = 1)
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
  ############

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
        if (SHOW == TRUE) {
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
#' @param SHOW SHOW
#' @param COL COL
#'
#' @return output
#' @export
vector.buildNet <- function(
    OUT,
    CUT = 1,
    SHOW = TRUE,
    COL = "grey70") {
  VEC <- OUT$VEC
  INDEX_LIST <- OUT$INDEX_LIST
  CENTER_LIST <- OUT$CENTER_LIST
  this_step_x <- OUT$this_step_x
  this_step_y <- OUT$this_step_y
  delta <- 0.00001

  if (SHOW == TRUE) {
    graphics::plot(
      VEC,
      col = COL,
      pch = 16,
      cex = 0.2,
      xlab = "",
      ylab = "",
      axes = FALSE,
      asp = 1
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
        if (SHOW == TRUE) {
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
  if (SHOW == TRUE) {
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
#' @param SHOW SHOW
#'
#' @return output
#' @export
vector.getValue <- function(OUT, PCA, SHOW = TRUE) {
  VALUE.OUT <- vector.calValue(PCA)

  OUT$VALUE <- VALUE.OUT$VALUE
  OUT$PCA <- PCA
  OUT$PCA.RC <- VALUE.OUT$PCA.RC

  if (SHOW == TRUE) {
    vector.showValue(OUT)
  }

  return(OUT)
}

#' Title
#'
#' @param OUT OUT
#' @param SHOW SHOW
#'
#' @return output
#' @export
vector.gridValue <- function(OUT, SHOW = TRUE) {
  INDEX_LIST <- OUT$INDEX_LIST
  VALUE <- OUT$VALUE
  SHOW <- SHOW
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

  if (SHOW == TRUE) {
    graphics::plot(OUT$VEC, col = "grey80", pch = 16, cex = 0.5, xlab = "", ylab = "", axes = FALSE, asp = 1)
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
#' @param SHOW SHOW
#'
#' @return output
#' @export
vector.autoCenter <- function(OUT, UP = 0.9, SHOW = TRUE) {
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

  if (SHOW == TRUE) {
    # plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
    graphics::plot(OUT$VEC, col = OUT$ORIG.COL, pch = 16, cex = 0.5, xlab = "", ylab = "", axes = FALSE, asp = 1)
    graphics::rect(
      graphics::par("usr")[1],
      graphics::par("usr")[3],
      graphics::par("usr")[2],
      graphics::par("usr")[4],
      col = NA,
      border = "black")
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

#' Title
#'
#' @param OUT OUT
#' @param P P
#' @param SHOW SHOW
#' @param COL COL
#' @param OL OL
#' @param AL AL
#' @param CEX CEX
#' @param AW AW
#' @param BD BD
#' @param AC AC
#' @param SHOW.SUMMIT SHOW.SUMMIT
#'
#' @return output
#' @export
vector.drawArrow <- function(
    OUT,
    P = 0.9,
    SHOW = TRUE,
    COL = "grey70",
    OL = 1.5,
    AL = 60,
    CEX = 0.5,
    AW = 2,
    BD = TRUE,
    AC = "grey20",
    SHOW.SUMMIT = TRUE) {
  USED <- OUT$USED
  DIST <- OUT$DIST
  ALL_VEC <- OUT$VEC
  USED_NAME <- OUT$USED_NAME
  USED_CENTER_VEC <- OUT$CENTER_VEC[USED, ]
  USED_DIST <- OUT$DIST[which(rownames(DIST) %in% USED_NAME), which(rownames(DIST) %in% USED_NAME)]
  SCORE <- OUT$SCORE

  one <- min(stats::dist(USED_CENTER_VEC)) * OL

  DIV <- 1 / P

  if (SHOW == TRUE) {
    # plot(ALL_VEC,col=COL,pch=16,cex=CEX)
    if (BD == TRUE) {
      graphics::plot(
        ALL_VEC,
        col = COL, pch = 16, cex = CEX,
        xlim = c(min(ALL_VEC[, 1]) - one, max(ALL_VEC[, 1]) + one),
        ylim = c(min(ALL_VEC[, 2]) - one, max(ALL_VEC[, 2]) + one),
        xlab = "", ylab = "", axes = FALSE, asp = 1
      )
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NA,
        border = "black")
    } else {
      graphics::plot(
        ALL_VEC,
        col = COL, pch = 16, cex = CEX,
        xlim = c(min(ALL_VEC[, 1]) - one, max(ALL_VEC[, 1]) + one),
        ylim = c(min(ALL_VEC[, 2]) - one, max(ALL_VEC[, 2]) + one), yaxt = "n", axes = F,
        xlab = "", ylab = "", axes = FALSE, asp = 1
      )
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NA,
        border = "black")
    }
  }
  N.SCORE <- .normX(SCORE)
  SCORE.COL <- vector.vcol(N.SCORE, c(0, 0.5, 1), c("#009FFF", "#FFF200", "#ec2F4B"))

  A1_VEC <- c()
  A2_VEC <- c()
  A_LENGTH <- c()
  # one=min(dist(USED_CENTER_VEC)) * OL

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

    # plot(ALL_VEC,col=COL,pch=16,cex=0.2)
    # this_arrow_length=0.1*sqrt(sum(final_vec^2))
    # this_arrow_length=sqrt(sum(final_vec^2))#0.1 #* (1+sqrt(sum(final_vec^2)))
    this_arrow_length <- grDevices::dev.size()[1] / AL * sqrt(sum(final_vec^2)) / one # * sqrt(sum(final_vec^2)) #0.25
    if (SHOW == TRUE) {
      this_arrow_col <- AC # grey20'
      # if(AC==TRUE){this_arrow_col=SCORE.COL[i]}
      graphics::arrows(
        x0 = this_p1_loc[1], y0 = this_p1_loc[2],
        x1 = this_p2_loc[1], y1 = this_p2_loc[2],
        lwd = AW, length = this_arrow_length,
        col = this_arrow_col
        # col=this_arrow_col#SCORE.COL[i]
      )
    }
    A1_VEC <- cbind(A1_VEC, this_p1_loc)
    A2_VEC <- cbind(A2_VEC, this_p2_loc)
    A_LENGTH <- c(A_LENGTH, this_arrow_length)
    i <- i + 1
  }

  if (SHOW == TRUE & SHOW.SUMMIT == TRUE) {
    # graphics::points(OUT$CENTER_VEC[OUT$SUMMIT,],col='black',pch=16,cex=1)
    # graphics::points(OUT$CENTER_VEC[OUT$SUMMIT,],col='red',pch=16,cex=0.5)
    X1 <- min(OUT$CENTER_VEC[OUT$SUMMIT, 1]) - one / 10
    X2 <- max(OUT$CENTER_VEC[OUT$SUMMIT, 1]) + one / 10
    Y1 <- min(OUT$CENTER_VEC[OUT$SUMMIT, 2]) - one / 10
    Y2 <- max(OUT$CENTER_VEC[OUT$SUMMIT, 2]) + one / 10

    graphics::rect(
      xleft = X1, ybottom = Y1, xright = X2, ytop = Y2, angle = 45,
      col = NA, border = "#009FFF", lty = 1, lwd = 2
    )
    # rect(xleft=X1, ybottom=Y1, xright=X2, ytop=Y2, angle = 45,
    # col = NA, border = 'black', lty = 2, lwd = 2)
  }
  A1_VEC <- t(A1_VEC)
  A2_VEC <- t(A2_VEC)
  # OUT=list()
  OUT$A1_VEC <- A1_VEC
  OUT$A2_VEC <- A2_VEC
  OUT$A_LENGTH <- A_LENGTH
  OUT$A_COL <- SCORE.COL

  return(OUT)
}
