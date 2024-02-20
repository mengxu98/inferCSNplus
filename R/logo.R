#' @title inferCSN logo
#'
#' @description
#' The inferCSN logo, using ASCII or Unicode characters
#' Use [cli::ansi_strip()] to get rid of the colors.
#' Reference: https://github.com/tidyverse/tidyverse/blob/main/R/logo.R
#' @param unicode Whether to use Unicode symbols. Default is `TRUE` on UTF-8 platforms.
#'
#' @md
#' @export
#' @examples
#' inferCSN.logo()
inferCSN.logo <- function(
    unicode = cli::is_utf8_output()) {
  logo <- c(
    "           2   1                3
            ____     4    ___________ _   __  __
    0 ___  / __/__  _____/ ____/ ___// | / /_/ /_
  / / __ ./ /_/ _ ./ ___/ /    |__ ./  |/ /_  __/
 / / / / / __/  __/ /  / /___ ___/ / /|  / /_/
/_/_/ /_/_/  .___/_/   .____//____/_/ |_/
      6      5               7      8       9   ")

  hexa <- c("*", ".", "o", "*", ".", "*", ".", "o", ".", "*")
  if (unicode) hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]

  cols <- c("red", "yellow", "green", "magenta", "cyan",
            "yellow", "green", "white", "magenta", "cyan")

  col_hexa <- purrr::map2(hexa, cols, ~ cli::make_ansi_style(.y)(.x))

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(cli::col_blue(logo), class = "logo")
}


#' print logo figure
#' @param x Input
#' @param ... Others
#'
#' @method print logo
#'
#' @export
print.logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}
