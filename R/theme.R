#' Removes legend from plot.
#'
#' @export
no_legend <- function() {
  theme(
    legend.position = "none"
  )
}

#' Removes margins from plot.
#'
#' @export
no_margin <- function() {
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "lines")
  )
}


#' Removes x axis text.
#'
#' @export
no_x_text <- function() {
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
}

#' Removes y axis text.
#'
#' @export
no_y_text <- function() {
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
}
