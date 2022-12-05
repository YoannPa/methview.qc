#' Checks if x is null and if so uses y by default.
#' 
#' @param x A value to be checked for being NULL or not.
#' @param y Another value to be returned if x is NULL.
#' @return A value that is not NULL.
#' @author Yoann Pageaud.
#' @references Henry L, Wickham H (2022). rlang: Functions for Base Types and
#'             Core R and 'Tidyverse' Features. https://rlang.r-lib.org,
#'             https://github.com/r-lib/rlang.
#' @keywords internal

null_default <- function(x, y){
  if (is.null(x)){ y } else { x }
}

#' Fixes display issues of geom_boxplot and geom_crossbar legend keys.
#' 
#' @author Yoann Pageaud.
#' @keywords internal

draw_key_boxplot2 <- function(data, params, size) {
  gp <- gpar(
    col = null_default(data$colour, "grey20"),
    fill = alpha(null_default(data$fill, "white"), data$alpha),
    lwd = null_default(data$linewidth, 0.5) * .pt,
    lty = null_default(data$linetype, 1),
    lineend = null_default(params$lineend, "butt"),
    linejoin = null_default(params$linejoin, "mitre")
  )
  if (isTRUE(params$flipped_aes)) {
    grobTree(
      linesGrob(c(0.1, 0.25), 0.5),
      linesGrob(c(0.75, 0.9), 0.5),
      rectGrob(width = 0.5, height = 0.75),
      linesGrob(0.5, c(0.125, 0.875)),
      gp = gp
    )
  } else {
    grobTree(
      linesGrob(0.5, c(0.1, 0.25)),
      linesGrob(0.5, c(0.75, 0.9)),
      rectGrob(height = 0.5, width = 0.75),
      linesGrob(c(0.125, 0.875), 0.5),
      gp = gp
    )
  }
}