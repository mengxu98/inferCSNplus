#' @title density_points
#'
#' @param meta_data meta_data
#' @param cluster_list cluster_list
#' @param method method
#' @param cluster_by cluster_by
#' @param group_column group_column
#' @param pseudotime_column pseudotime_column
#' @param min_cells min_cells
#' @param plot plot
#' @param key_word key_word
#' @param key_sign key_sign
#' @param file_save file_save
#' @param color_list color_list
#'
#' @return points
#' @export
#'
#' @examples
#' test_data <- rbind(
#'   data.frame(
#'     cluster = "cluster1",
#'     pseudotime = rnorm(500, mean = 1, sd = 1)
#'   ),
#'   data.frame(
#'     cluster = "cluster2",
#'     pseudotime = rnorm(500, mean = 2, sd = 2)
#'   ),
#'   data.frame(
#'     cluster = "cluster3",
#'     pseudotime = rnorm(500, mean = 3, sd = 3)
#'   )
#' )
#'
#' density_points(
#'   meta_data = test_data,
#'   plot = TRUE
#' )
#'
#' density_points(
#'   meta_data = test_data,
#'   cluster_list = list(c("cluster1", "cluster3")),
#'   plot = TRUE,
#'   color_list = c("#ff0000", "#ffcc00")
#' )
#'
#' density_points(
#'   meta_data = test_data,
#'   method = "intersection",
#'   plot = TRUE
#' )
density_points <- function(
    meta_data,
    cluster_list = NULL,
    method = "combine",
    cluster_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime",
    min_cells = 100,
    plot = FALSE,
    key_word = "network_",
    key_sign = "..",
    file_save = NULL,
    color_list = c(
      "#0066ff", "#0099ff", "#66cc33", "#66ff33", "#ccff00", "#ffff99",
      "#ffff33", "#ffcc33", "#ff9933", "#ff6633", "#cc3300", "#ff3333"
    )) {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]

  meta_data$cluster <- meta_data[, group_column]

  meta_data <- meta_data[, c("cluster", "pseudotime")]

  if (is.null(cluster_list)) {
    cluster_by <- unique(meta_data$cluster)
  } else{
    cluster_by <- unique(unlist(cluster_list))
  }
  meta_data <- purrr::map_dfr(
    cluster_by, function(x) {
      res <- dplyr::filter(meta_data, cluster == x)
      if (nrow(res) >= min_cells) {
        return(res)
      } else {
        return()
      }
    }
  )

  res <- purrr::map_dfr(
    unique(meta_data$cluster), function(x) {
      meta_data_sub <- dplyr::filter(meta_data, cluster == x)
      pseudotime <- meta_data_sub$pseudotime
      xlim <- c(min(pseudotime), max(pseudotime))
      df <- as.data.frame(
        stats::density(
          pseudotime,
          kernel = "gaussian",
          n = length(pseudotime),
          from = xlim[1],
          to = xlim[2]
        )[c("x", "y")]
      )
      max_y <- max(df$y)
      x_point <- df$x[which(df$y == max(df$y))]
      data.frame(
        cluster = x,
        max_density = max_y,
        max_density_pseudotime = x_point
      )
    }
  )
  res <- res[order(res$max_density_pseudotime), ]
  cluster_order <- res$cluster

  if (is.null(cluster_list)) {
    cluster_list <- list()
    for (i in seq_len(nrow(res) - 1)) {
      cluster_list[[i]] <- c(res$cluster[i], res$cluster[i + 1])
    }
  }

  if (method == "combine") {
    points_data <- res

    intersection_points <- as.numeric(points_data$max_density_pseudotime)
    all_intersection_points <- intersection_points

    results <- data.frame()
    for (i in seq_len(nrow(points_data) - 1)) {
      results <- rbind(
        results,
        data.frame(
          window = paste0(key_word, i),
          cluster1 = points_data$cluster[i],
          cluster2 = points_data$cluster[i + 1],
          time_point1 = points_data$max_density_pseudotime[i],
          time_point2 = points_data$max_density_pseudotime[i + 1]
        )
      )
    }
  }

  if (method == "intersection") {
    points_data <- purrr::map_dfr(
      cluster_list, function(x) {
        a <- dplyr::filter(meta_data, cluster == x[1])$pseudotime
        b <- dplyr::filter(meta_data, cluster == x[2])$pseudotime

        xlim <- c(min(c(a, b)), max(c(a, b)))
        df <- merge(
          as.data.frame(
            stats::density(
              a,
              kernel = "gaussian",
              from = xlim[1],
              to = xlim[2]
            )[c("x", "y")]
          ),
          as.data.frame(
            stats::density(
              b,
              kernel = "gaussian",
              from = xlim[1],
              to = xlim[2]
            )[c("x", "y")]
          ),
          by = "x",
          suffixes = c(".a", ".b")
        )

        df$comp <- as.numeric(df$y.a > df$y.b)
        df$cross <- c(NA, diff(df$comp))
        intersection_point <- df[which(df$cross != 0), "x"]
        if (length(intersection_point) > 1) {
          intersection_point <- df[which(df$cross == (-1)), "x"]
          intersection_point <- max(intersection_point)
        }
        data.frame(
          cluster1 = x[1],
          cluster2 = x[2],
          max_density = max(c(max(df$y.a), max(df$y.b))),
          max_density_pseudotime = intersection_point
        )
      }
    )

    intersection_points <- as.numeric(points_data$max_density_pseudotime)
    all_intersection_points <- c(
      min(meta_data$pseudotime),
      intersection_points,
      max(meta_data$pseudotime)
    )

    points_data <- rbind(
      data.frame(
        cluster1 = cluster_list[[1]][1],
        cluster2 = cluster_list[[1]][1],
        max_density = 0,
        max_density_pseudotime = min(meta_data$pseudotime)
      ),
      points_data,
      data.frame(
        cluster1 = cluster_list[[length(cluster_list)]][2],
        cluster2 = cluster_list[[length(cluster_list)]][2],
        max_density = 0,
        max_density_pseudotime = max(meta_data$pseudotime)
      )
    )

    results <- data.frame()
    for (i in seq_len(nrow(points_data) - 2)) {
      results <- rbind(
        results,
        data.frame(
          window = paste0(key_word, i),
          cluster1 = points_data$cluster1[i],
          cluster2 = points_data$cluster2[i + 2],
          time_point1 = points_data$max_density_pseudotime[i],
          time_point2 = points_data$max_density_pseudotime[i + 2]
        )
      )
    }
  }

  if (plot) {
    max_density_point <- max(points_data$max_density)

    plot_data <- meta_data[, c("cluster", "pseudotime")]
    plot_data$cluster <- factor(
      plot_data$cluster,
      levels = cluster_order,
      labels = cluster_order
    )

    xIndex <- split(plot_data[["pseudotime"]], plot_data[["cluster"]])
    bg_data <- as.data.frame(t(sapply(xIndex, range)))
    colnames(bg_data) <- c("xmin", "xmax")
    bg_data[["group.by"]] <- names(xIndex)

    bg_data[["ymin"]] <- 0
    bg_data[["ymax"]] <- Inf
    bg_data[["fill"]] <- color_list[seq_along(cluster_order)]

    xmin <- c()
    xmax <- c()
    if (length(cluster_order) > 1) {
      for (x in seq_along(cluster_order)) {
        if (x == 1) {
          xmin[x] <- min(plot_data$pseudotime)
          xmax[x] <- intersection_points[x]
        } else if (x == length(cluster_order)) {
          xmin[x] <- intersection_points[x - 1]
          xmax[x] <- max(plot_data$pseudotime)
        } else {
          xmin[x] <- intersection_points[x - 1]
          xmax[x] <- intersection_points[x]
        }
      }
    }
    bg_data[["xmin"]] <- xmin
    bg_data[["xmax"]] <- xmax

    bg_layer <- geom_rect(
      data = bg_data,
      xmin = bg_data[["xmin"]], xmax = bg_data[["xmax"]],
      ymin = bg_data[["ymin"]], ymax = bg_data[["ymax"]],
      fill = bg_data[["fill"]],
      alpha = 0.2,
      inherit.aes = FALSE
    )

    p <- ggplot(data = plot_data, aes(x = pseudotime))
    p <- p + bg_layer
    p <- p +
      geom_density(aes(color = cluster)) +
      geom_density(aes(fill = cluster), alpha = 0.8) +
      geom_vline(
        xintercept = intersection_points,
        color = "#df1d79",
        linetype = "dashed"
      ) +
      scale_fill_manual(values = color_list) +
      labs(x = "Pseudotime", y = "Density") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_bw()
    p <- p + theme(legend.position = "bottom")

    if (method == "intersection") {
      y_points <- c()
      x_window <- c()
      y_window <- c()
      for (w in seq_along(cluster_list)) {
        y_points <- c(y_points, max_density_point * w / length(cluster_order))
        x_window <- c(
          x_window,
          mean(c(
            all_intersection_points[w],
            all_intersection_points[w + 2]
          ))
        )
        y_window <- c(y_window, max_density_point * w / length(cluster_order) + 0.1)
      }

      line_data <- data.frame(
        x = all_intersection_points[seq_along(cluster_list)],
        xend = all_intersection_points[(length(all_intersection_points) - length(cluster_list) + 1):length(all_intersection_points)],
        y = sort(y_points, decreasing = TRUE),
        x_window = x_window,
        y_window = sort(y_window, decreasing = TRUE),
        window = paste0(key_word, seq_along(cluster_list))
      )
    }

    if (method == "combine") {
      y_points <- c()
      x_window <- c()
      y_window <- c()
      for (w in seq_along(cluster_list)) {
        y_points <- c(y_points, max_density_point * w / length(cluster_order))
        x_window <- c(
          x_window,
          mean(c(
            all_intersection_points[w],
            all_intersection_points[w + 1]
          ))
        )
        y_window <- c(y_window, max_density_point * w / length(cluster_order) + 0.1)
      }

      line_data <- data.frame(
        x = all_intersection_points[seq_along(cluster_list)],
        xend = all_intersection_points[(length(all_intersection_points) - length(cluster_list) + 1):length(all_intersection_points)],
        y = sort(y_points, decreasing = TRUE),
        x_window = x_window,
        y_window = sort(y_window, decreasing = TRUE),
        window = paste0(key_word, seq_along(cluster_list), key_sign, (seq_along(cluster_list) + 1))
      )
    }

    p <- p + geom_segment(
      data = line_data,
      aes(
        x = x, xend = xend,
        y = y, yend = y
      ),
      color = "gray70",
      linewidth = 1.3,
      arrow = arrow(angle = 15, type = "closed")
    ) +
      geom_text(
        data = line_data,
        # linewidth = 5,
        aes(
          x = x_window,
          y = y_window,
          label = window
        )
      )

    print(p)
    if (!is.null(file_save)) {
      cowplot::ggsave2(
        file = file_save,
        p,
        width = 15,
        height = 6,
        units = "cm",
        dpi = 600
      )
    }
  }

  return(results)
}

#' dynamic.windowing
#'
#' @param meta_data meta_data
#' @param cluster_list cluster_list
#' @param bin_points bin_points
#' @param cluster_by cluster_by
#' @param group_column group_column
#' @param pseudotime_column pseudotime_column
#' @param min_cells min_cells
#' @param key_word key_word
#' @param key_sign key_sign
#'
#' @return a
#' @export
dynamic.windowing <- function(
    meta_data,
    bin_points = NULL,
    cluster_list = NULL,
    cluster_by = NULL,
    group_column = "cluster",
    pseudotime_column = "pseudotime",
    min_cells = 100,
    key_word = "network_",
    key_sign = "..") {
  meta_data$pseudotime <- meta_data[, pseudotime_column]
  meta_data <- meta_data[which(!is.na(meta_data$pseudotime)), ]

  meta_data$cluster <- meta_data[, group_column]
  meta_data$cells <- rownames(meta_data)

  meta_data <- meta_data[, c("cluster", "pseudotime", "cells")]

  if (is.null(cluster_list)) {
    cluster_by <- unique(meta_data$cluster)
  } else{
    cluster_by <- unique(unlist(cluster_list))
  }
  meta_data <- purrr::map_dfr(
    cluster_by, function(x) {
      res <- filter(meta_data, cluster == x)
      if (nrow(res) >= min_cells) {
        return(res)
      }
    }
  )
  cluster_by <- unique(meta_data$cluster)

  res <- purrr::map_dfr(
    cluster_by, function(x) {
      meta_data_sub <- dplyr::filter(meta_data, cluster == x)
      pseudotime <- meta_data_sub$pseudotime
      xlim <- c(min(pseudotime), max(pseudotime))
      df <- as.data.frame(
        stats::density(
          pseudotime,
          kernel = "gaussian",
          n = length(pseudotime),
          from = xlim[1],
          to = xlim[2]
        )[c("x", "y")]
      )
      max_y <- max(df$y)
      x_point <- df$x[which(df$y == max(df$y))]
      data.frame(
        cluster = x,
        max_density = max_y,
        max_density_pseudotime = x_point
      )
    }
  )
  res <- res[order(res$max_density_pseudotime), ]
  list_names1 <- c()
  for (i in seq_len(nrow(res))) {
    list_names1 <- c(list_names1, paste0(key_word, i, key_sign, i))
  }

  cells_list1 <- apply(
    res, 1, function(x) {
      dplyr::filter(meta_data, cluster == x[1])
    }
  )
  names(cells_list1) <- list_names1

  if (is.null(bin_points)) {
    bin_points <- density_points(
      meta_data,
      cluster_list = cluster_list
    )
  }

  list_names2 <- c()
  for (i in seq_len(nrow(bin_points))) {
    list_names2 <- c(list_names2, paste0(key_word, i, key_sign, i + 1))
  }

  cells_list2 <- apply(
    bin_points, 1, function(x) {
      meta_data_sub1 <- dplyr::filter(meta_data, cluster == x[2])
      meta_data_sub2 <- dplyr::filter(meta_data, cluster == x[3])
      meta_data_sub1 <- meta_data_sub1[which(meta_data_sub1$pseudotime > x[4]), ]
      meta_data_sub2 <- meta_data_sub2[which(meta_data_sub2$pseudotime <= x[5]), ]
      rbind(
        meta_data_sub1,
        meta_data_sub2
      )
    }
  )
  names(cells_list2) <- list_names2

  cells_list <- c(cells_list1, cells_list2)

  return(cells_list)
}
