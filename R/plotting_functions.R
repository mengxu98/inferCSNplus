#' @title plot_scatter
#'
#' @param data data
#' @param method method
#' @param group_colors group_colors
#' @param titleColor titleColor
#' @param title title
#' @param xTitle xTitle
#' @param yTitle yTitle
#' @param legendTitle legendTitle
#' @param legend legend
#'
#' @return ggplot object
#' @export
#' @examples
#' data("example_matrix")
#' test_data <- data.frame(
#'   example_matrix[, c(1,11)],
#'   c = c(rep("c1", 2000), rep("c2", 2000), rep("c3", 1000))
#' )
#' plot_scatter(test_data)
plot_scatter <- function(
    data,
    method = "loess",
    group_colors = c("#e1181f", "#255830", "#006699"),
    titleColor = "#006699",
    title = NULL,
    xTitle = NULL,
    yTitle = NULL,
    legendTitle = NULL,
    legend = "bottom") {
  method <- match.arg(method, c("lm", "loess"))
  if (ncol(data) == 3) {
    colnames(data) <- c("x", "y", "cluster")
    p <- ggplot(data = data, mapping = aes(x = x, y = y, cluster = cluster, col = cluster)) +
      geom_smooth(formula = 'y ~ x',
                  method = method)
  } else if (ncol(data) == 2) {
    colnames(data) <- c("x", "y")
    group_colors <- group_colors[1]
    p <- ggplot(data = data, mapping = aes(x = x, y = y), col = "gray") +
      geom_smooth(formula = 'y ~ x',
                  color = group_colors,
                  method = method)
  }

  p <- p +
    geom_point() +
    theme_bw() +
    ggpubr::stat_cor(data = data, method = "pearson") +
    labs(title = title, x = xTitle, y = yTitle, fill = legendTitle) +
    theme(plot.title = element_text(color = titleColor),
          legend.position = legend) +
    scale_color_manual(values = group_colors)

  # ggpubr::ggscatter(data,
  #                   "x", "y",
  #                   color = "#00AFBB",
  #                   # size = "qsec",
  #                   add = "loess",
  #                   conf.int = TRUE) + theme_bw()
  # smoothScatter(data$x, data$y)
  return(p)
}
