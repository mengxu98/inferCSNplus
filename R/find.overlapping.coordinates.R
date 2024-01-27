#' peaks.split
#'
#' @description
#' This function from:\url{https://github.com/cole-trapnell-lab/cicero-release/blob/master/R/utils.R#L507}
#'
#' @param inp parameter
#'
#' @return seq
#' @export
peaks.split <- function(inp) {
  out <- stringr::str_split_fixed(
    stringi::stri_reverse(inp),
    ":|-|_", 3
  )
  out[, 1] <- stringi::stri_reverse(out[, 1])
  out[, 2] <- stringi::stri_reverse(out[, 2])
  out[, 3] <- stringi::stri_reverse(out[, 3])
  out[, c(3, 2, 1), drop = FALSE]
}

#' Construct GRanges objects from coordinate strings
#'
#' @description
#' This function from:\url{https://github.com/cole-trapnell-lab/cicero-release/blob/master/R/utils.R#L90}
#'
#' @param coord_strings A list of coordinate strings (in the form "chr1:500000-1000000").
#' @param with_names logical - should meta data include coordinate string (field coord_string)?
#' @param meta_data_df A data frame with any meta data columns you want
#'   included with the ranges. Must be in the same order as coord_strings.
#'
#' @details Coordinate strings consist of three pieces of information:
#'   chromosome, start, and stop. These pieces of information can be separated
#'   by the characters ":", "_", or "-". Commas will be removed, not used as
#'   separators (ex: "chr18:8,575,097-8,839,855" is ok).
#'
#' @return GRanges object of the input strings
#' @export
#'
#' @seealso \code{\link[GenomicRanges]{GRanges-class}}
#'
#' @examples
#' ran1 <- ranges.for.coords("chr1:2039-30239", with_names = TRUE)
#' ran2 <- ranges.for.coords(c("chr1:2049-203902", "chrX:489249-1389389"),
#'   meta_data_df = data.frame(dat = c("1", "X"))
#' )
#' ran3 <- ranges.for.coords(c("chr1:2049-203902", "chrX:489249-1389389"),
#'   with_names = TRUE,
#'   meta_data_df = data.frame(
#'     dat = c("1", "X"),
#'     stringsAsFactors = FALSE
#'   )
#' )
ranges.for.coords <- function(
    coord_strings,
    meta_data_df = NULL,
    with_names = FALSE) {
  assertthat::assert_that(is.logical(with_names))
  if (!is.null(meta_data_df)) {
    assertthat::assert_that(is.data.frame(meta_data_df))
    assertthat::assert_that(
      assertthat::are_equal(length(coord_strings), nrow(meta_data_df))
    )
  }

  coord_strings <- gsub(",", "", coord_strings)
  coord_cols <- peaks.split(coord_strings)
  gr <- GenomicRanges::GRanges(
    coord_cols[, 1],
    ranges = IRanges::IRanges(as.numeric(coord_cols[, 2]), as.numeric(coord_cols[, 3])),
    mcols = meta_data_df
  )
  if (!is.null(meta_data_df)) {
    for (n in names(meta_data_df)) {
      newname <- paste0("mcols.", n)
      names(GenomicRanges::mcols(gr))[which(names(GenomicRanges::mcols(gr)) == newname)] <- n
    }
  }
  if (with_names) {
    gr$coord_string <- coord_strings
  }
  gr
}

#' Find peaks that overlap a specific genomic location
#'
#' @description
#' This function from:\url{https://github.com/cole-trapnell-lab/cicero-release/blob/master/R/utils.R#L438}
#'
#' @param coord_list A list of coordinates to be searched for overlap in the
#'   form chr_100_2000.
#' @param coord The coordinates that you want to find in the form chr1_100_2000.
#' @param maxgap The maximum distance in base pairs between coord and the
#'  coord_list that should count as overlapping. Default is 0.
#'
#' @return A character vector of the peaks that overlap coord.
#' @export
#'
#' @examples
#' test_coords <- c(
#'   "chr18_10025_10225", "chr18_10603_11103",
#'   "chr18_11604_13986",
#'   "chr18_157883_158536", "chr18_217477_218555",
#'   "chr18_245734_246234"
#' )
#' find.overlapping.coordinates(test_coords, "chr18:10,100-1246234")
find.overlapping.coordinates <- function(
    coord_list,
    coord,
    maxgap = 0) {
  coord <- gsub(",", "", coord)
  cons_gr <- ranges.for.coords(coord_list)
  if (length(coord) == 1) {
    ol1 <- GenomicRanges::findOverlaps(
      ranges.for.coords(coord),
      cons_gr,
      maxgap = maxgap,
      select = "all"
    )
    ol1 <- as.list(ol1)
    return(as.character(coord_list[unlist(ol1)]))
  } else {
    ol1 <- lapply(coord, function(x) {
      y <- suppressWarnings(
        unlist(as.list(
          GenomicRanges::findOverlaps(ranges.for.coords(x),
            cons_gr,
            maxgap = maxgap,
            select = "all"
          )
        ))
      )
      if (length(y) == 0) {
        return(NA)
      }
      return(coord_list[y])
    })

    return(as.character(unlist(ol1)))
  }
}
