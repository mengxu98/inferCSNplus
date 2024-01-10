#' extract links
#' @param cre_gene_list A list of CRE-gene
#' @param cre_tf_list A list of CRE-TF
#' @param promoter_cre_gene_list A list of CRE-gene about promoters
#' @param promoter_cre_tf_list A list of CRE-TF about promoters
#' @param clusters cluster information
#'
#' @return TF-G dataframe
#' @export
extract.links <- function(
    cre_gene_list,
    cre_tf_list,
    promoter_cre_gene_list,
    promoter_cre_tf_list,
    clusters) {
  tf_gene_cre <- get.tf.gene(cre_gene_list, cre_tf_list, clusters)
  promoter_tf_gene_cre <- get.tf.gene(promoter_cre_gene_list, promoter_cre_tf_list, clusters)
  # take genes which has both distal TFs and proximal TFs
  intersect_genes <- intersect(tf_gene_cre$Gene, promoter_tf_gene_cre$Gene)
  tf_gene_cre <- tf_gene_cre[which(tf_gene_cre$Gene %in% intersect_genes), , drop = FALSE]
  promoter_tf_gene_cre <- promoter_tf_gene_cre[which(promoter_tf_gene_cre$Gene %in% intersect_genes), , drop = FALSE]
  tf_gene_cre_final <- rbind(tf_gene_cre, promoter_tf_gene_cre)
  tf_gene_cre_final$types <- c(
    rep("distal", nrow(tf_gene_cre)),
    rep("proximal", nrow(promoter_tf_gene_cre))
  )
  return(tf_gene_cre_final)
}

#' extract Node character
#' @param network_links Links obtained by extract.links
#' @param markers Dataframe with marker and cluster information
#' @return Node character dataframe
#' @export
extract.nodes <- function(
    network_links,
    markers) {
  ### Node
  G <- unique(network_links$Gene)
  T <- unique(network_links$TF)
  Nodes <- c(T, G)
  label_type <- c(rep("TF", length(T)), rep("gene", length(G)))
  unik <- !duplicated(Nodes)
  Nodes <- Nodes[unik]
  label_type <- label_type[unik]

  colnames(markers) <- c("gene", "cluster")
  names <- unique(markers$cluster)

  ### used to note cluster specific information
  DE_flag <- matrix(0, nrow = length(Nodes), ncol = length(names))
  for (i in 1:length(names)) {
    marker1 <- markers[markers$cluster == names[i], ]
    marker_gene <- as.character(marker1$gene)
    id <- which(Nodes %in% marker_gene)
    DE_flag[id, i] <- 1
  }
  cluster_type <- rep(NA, length(Nodes))
  for (i in 1:length(Nodes)) {
    id <- which(DE_flag[i, ] == 1)
    name_i <- names[id[1]]
    if (length(id) == 1) {
      cluster_type[i] <- name_i
    } else if (length(id) == 2) {
      cluster_type[i] <- paste0(name_i, "_", names[id[2]])
    } else {
      for (j in 2:(length(id) - 1)) {
        name_i <- paste0(name_i, "_", names[id[j]])
      }
      cluster_type[i] <- paste0(name_i, "_", names[id[length(id)]])
    }
  }
  Node_character <- cbind(label_type, cluster_type)
  colnames(Node_character) <- c("function_type", "cluster")
  rownames(Node_character) <- Nodes
  return(Node_character)
}

#' Obtain TF and Gene relationships
#' @param cre_gene_list A list of CRE-gene
#' @param cre_tf_list A list of CRE-TF
#' @param clusters cluster information
#' @return TF-G dataframe
get.tf.gene <- function(
    cre_gene_list,
    cre_tf_list,
    clusters) {
  L_TF_G_record <- list()
  tf_gene_cre_record <- list()
  for (i in 1:length(cre_tf_list)) {
    L_T <- cre_tf_list[[i]]
    L_G <- cre_gene_list[[i]]
    genes <- rep(NA, nrow(L_T))
    for (j in 1:nrow(L_T)) {
      id <- which(L_G$loci %in% L_T$loci[j])
      if (length(id) > 0) {
        genes[j] <- L_G$gene[id]
      }
    }
    L_TF_G <- L_T[!is.na(genes), , drop = FALSE]
    L_TF_G$gene <- genes[!is.na(genes)]
    L_TF_G_record[[i]] <- L_TF_G

    TF_G <- L_TF_G[, 2:3]
    rownames(TF_G) <- NULL
    # head(TF_G)
    unik <- !duplicated(TF_G)
    TF_G <- TF_G[unik, ]
    tf_gene_cre_record[[i]] <- data.frame(
      TF = TF_G$TF,
      Gene = TF_G$gene,
      Cell_type = clusters[i]
    )
  }
  tf_gene_cre <- do.call(rbind, tf_gene_cre_record)
  unik <- !duplicated(tf_gene_cre)
  tf_gene_cre <- tf_gene_cre[unik, ]
}

#' Extract cre-gene relations of markers
#'
#' @param weight_table The result of \code{inferCSN()}
#' @param markers A dataframe data with marker genes and cluster information,
#' such as the results of \code{Seurat::FindAllMarkers(object)}.
#' @param filter filter
#'
#' @export
extract.cre.gene.links <- function(
    weight_table,
    markers,
    filter = FALSE) {
  # loci_gene corresponding relationship
  if (filter) {
    weight_table <- weight_table[which(weight_table$function_type == "high_corr"), , drop = FALSE]
  }

  cre_gene_list <- list()
  promoter_cre_gene_list <- list() # promoter-gene
  uniq_clusters <- as.character(unique(markers$cluster))
  for (i in 1:length(uniq_clusters)) {
    marker1 <- markers[markers$cluster == uniq_clusters[i], ]
    marker_gene <- as.character(marker1$gene)
    effctive_id <- which(weight_table$gene %in% marker_gene)
    L_G_i <- data.frame(loci = weight_table$target[effctive_id], gene = weight_table$gene[effctive_id])
    cre_gene_list[[i]] <- L_G_i[!duplicated(L_G_i), , drop = FALSE]
    P_L_G_i <- data.frame(loci = weight_table$regulator[effctive_id], gene = weight_table$gene[effctive_id])
    promoter_cre_gene_list[[i]] <- P_L_G_i[!duplicated(P_L_G_i), , drop = FALSE]
  }
  CRE_Gene <- list()
  CRE_Gene$distal <- cre_gene_list
  CRE_Gene$promoter <- promoter_cre_gene_list
  return(CRE_Gene)
}

#' extract cre of markers
#' @param cre_gene_list a list of CRE-Gene relationships
#' @param promoter_cre_gene_list a list of Promoter-Gene relationships
#' @param da_peaks_list a list of DA of each cluster
#'
#' @export
extract.cre <- function(
    cre_gene_list,
    promoter_cre_gene_list,
    da_peaks_list = NULL) {
  # extract overlapped peaks between DA and cre of focused markers
  peaks_bed_record <- list() # peaks used to identify TFs bounded to
  P_peaks_bed_record <- list()

  if (is.null(da_peaks_list)) {
    for (i in 1:length(cre_gene_list)) {
      da_peaks_list[[i]] <- cre_gene_list[[i]]$loci
    }
  }
  cre_gene_list_new <- list()
  promoter_cre_gene_list_new <- list()
  for (i in 1:length(cre_gene_list)) {
    if (!is.null(da_peaks_list[[i]])) {
      if (!is.character(da_peaks_list[[i]])) {
        da <- rownames(da_peaks_list[[i]])
      } else {
        da <- da_peaks_list[[i]]
      }
      da <- gsub("-", "_", da)
    }

    DA_HC_i <- intersect(da, cre_gene_list[[i]]$loci)
    if (length(DA_HC_i) > 0) {
      cre_gene_list_new[[i]] <- cre_gene_list[[i]][which(cre_gene_list[[i]]$loci %in% DA_HC_i), ]
      overlap_gene <- cre_gene_list[[i]]$gene[which(cre_gene_list[[i]]$loci %in% DA_HC_i)]
      peaks1 <- strsplit(DA_HC_i, "_")
      peaks_bed <- matrix(0, nrow = length(DA_HC_i), ncol = 3)
      for (j in 1:length(DA_HC_i)) {
        for (k in 1:3) {
          peaks_bed[j, k] <- peaks1[[j]][k]
        }
      }
      colnames(peaks_bed) <- c("R.chrom", "R.start", "R.end")
      peaks_bed <- as.data.frame(peaks_bed)
      peaks_bed$R.start <- as.numeric(peaks_bed$R.start)
      peaks_bed$R.end <- as.numeric(peaks_bed$R.end)
      peaks_bed_record[[i]] <- peaks_bed
      # promoters
      promoter_cre_gene_list_i <- promoter_cre_gene_list[[i]]
      promoter_cre_gene_list_i <- promoter_cre_gene_list_i[which(promoter_cre_gene_list_i$gene %in% overlap_gene), ]
      promoter_cre_gene_list_new[[i]] <- promoter_cre_gene_list_i
      peaks2 <- strsplit(promoter_cre_gene_list_i$loci, "_")
      peaks_bed <- matrix(0, nrow = length(promoter_cre_gene_list_i$loci), ncol = 3)
      for (j in 1:length(promoter_cre_gene_list_i$loci)) {
        for (k in 1:3) {
          peaks_bed[j, k] <- peaks2[[j]][k]
        }
      }
      colnames(peaks_bed) <- c("R.chrom", "R.start", "R.end")
      peaks_bed <- as.data.frame(peaks_bed)
      peaks_bed$R.start <- as.numeric(peaks_bed$R.start)
      peaks_bed$R.end <- as.numeric(peaks_bed$R.end)
      P_peaks_bed_record[[i]] <- peaks_bed[!(duplicated(peaks_bed)), , drop = FALSE]
    }
  }
  cre <- list()
  cre$distal <- peaks_bed_record
  cre$promoter <- P_peaks_bed_record
  cre$cre_gene_list <- cre_gene_list_new
  cre$promoter_cre_gene_list <- promoter_cre_gene_list_new
  return(cre)
}

#' @title Identify TFs enriched in cre of focus markers
#'
#' @param peaks_bed_list A list of peaks bed file of each cluster
#' @param species Species used to detect TF
#' @param genome For example: BSgenome.Hsapiens.UCSC.hg38
#' @param markers two column dataframe data with marker genes and cluster information
#'
#' @export
extract.peak.tf.links <- function(
    peaks_bed_list,
    species,
    genome,
    markers) {
  motifs <- chromVAR::getJasparMotifs(species)
  cre_tf_list <- list()
  for (i in 1:length(peaks_bed_list)) {
    if (!is.null(peaks_bed_list[[i]])) {
      peaks_bed <- peaks_bed_list[[i]]
      peaks_new <- GenomicRanges::GRanges(
        seqnames = peaks_bed$R.chrom,
        ranges = IRanges::IRanges(
          start = peaks_bed$R.start,
          end = peaks_bed$R.end
        )
      )
      motif_ix <- motifmatchr::matchMotifs(
        motifs,
        peaks_new,
        genome,
        out = "scores"
      )
      S <- as.matrix(motif_ix@assays@data$motifScores)
      M <- as.matrix(motif_ix@assays@data$motifMatches)
      TF <- motif_ix@colData$name

      L_TF_list <- list()
      for (j in 1:nrow(M)) {
        if (sum(M[j, ]) > 0) {
          p <- paste0(peaks_bed$R.chrom[j], "_", peaks_bed$R.start[j], "_", peaks_bed$R.end[j])
          # focus TFs in the markers
          TF_j <- intersect(unique(markers$gene), TF[M[j, ]])
          if (length(TF_j) > 0) {
            L_TF_list[[j]] <- data.frame(loci = p, TF = TF_j)
          }
        }
      }
      cre_tf_list[[i]] <- do.call(rbind, L_TF_list)
    }
  }
  return(cre_tf_list)
}

#' @title extract.network
#'
#' @param object Seurat object with network results after run 'inferCSN()'
#' @param cluster Selected celltype(s)
#'
#' @return network
#' @export
extract.network <- function(
    object,
    cluster) {
  if ("weight_table_atac" %in% names(Seurat::Misc(object))) {
    weight_table <- Seurat::Misc(object, slot = "weight_table_atac")
  }
  weight_table <- as.data.frame(weight_table)

  if ("all_markers_list" %in% names(Seurat::Misc(object))) {
    all_markers_list <- Seurat::Misc(object, slot = "all_markers_list")
  }
  all_markers_list <- as.data.frame(all_markers_list)

  focused_markers <- all_markers_list[which(all_markers_list$cluster %in% cluster), , drop = FALSE]

  # CRE-gene connections
  cre_gene <- extract.cre.gene.links(
    weight_table,
    markers = focused_markers
  )

  # Find focused cre which is overlapped with DA
  cre <- extract.cre(
    cre_gene_list = cre_gene$distal,
    promoter_cre_gene_list = cre_gene$promoter
  )

  # detect TFs for distal cre
  cre_tf_list <- extract.peak.tf.links(
    peaks_bed_list = cre$distal,
    species = "Homo sapiens",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    markers = focused_markers
  )

  # detect TFs for Promoters
  promoter_cre_tf_list <- extract.peak.tf.links(
    peaks_bed_list = cre$promoter,
    species = "Homo sapiens",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    markers = focused_markers
  )

  network_links <- extract.links(
    cre_gene_list = cre_gene$distal,
    cre_tf_list = cre_tf_list,
    promoter_cre_gene_list = cre_gene$promoter,
    promoter_cre_tf_list = promoter_cre_tf_list,
    cluster
  )

  return(network_links)
}
