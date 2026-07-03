#' Fire-size distribution histogram with median-size-per-bin overlay
#'
#' Draws the fire-size distribution the way the `burnSummaries` module does: a
#' histogram of `log(fire size)` (the number of fires per size bin, left axis)
#' overlaid with the median `log(fire size)` within each bin (right/secondary
#' axis). Non-positive and non-finite sizes are dropped.
#'
#' @param sizes A `data.frame` of fire records (with column `size_col`) or a bare
#'   numeric vector of fire sizes. `NULL`, empty, or all-non-positive input returns
#'   `NULL`.
#' @param size_col Name of the fire-size column when `sizes` is a `data.frame`
#'   (default `"size"`).
#' @param size_unit Unit label for the size axis (default `"ha"`).
#' @param binwidth Width of the `log(size)` bins (default `0.5`).
#' @param title Optional plot title.
#'
#' @return A `ggplot` object, or `NULL` for empty input.
#'
#' @export
fire_size_histogram <- function(
  sizes,
  size_col = "size",
  size_unit = "ha",
  binwidth = 0.5,
  title = NULL
) {
  s <- if (is.data.frame(sizes)) sizes[[size_col]] else sizes
  s <- s[is.finite(s) & s > 0]
  if (!length(s)) {
    return(invisible(NULL))
  }
  logsize <- log(s)
  upper <- ceiling(max(logsize) / binwidth) * binwidth
  breaks <- seq(0, max(upper, binwidth), binwidth)
  bins <- cut(logsize, breaks, include.lowest = TRUE)
  counts <- as.integer(table(bins))
  med_log <- as.numeric(tapply(logsize, bins, stats::median))
  mids <- breaks[-length(breaks)] + binwidth / 2

  ## scale the median trace onto the count axis, then invert on the secondary axis
  max_count <- max(counts)
  max_med <- max(med_log, na.rm = TRUE)
  scale <- if (is.finite(max_med) && max_med > 0) max_count / max_med else 1

  bars <- data.frame(mid = mids, count = counts)
  pts <- data.frame(mid = mids, med = med_log)
  pts <- pts[is.finite(pts$med), , drop = FALSE]

  y2col <- "darkred"
  ggplot2::ggplot(bars, ggplot2::aes(x = .data[["mid"]], y = .data[["count"]])) +
    ggplot2::geom_col(width = binwidth, alpha = 0.5, fill = "grey20") +
    ggplot2::geom_point(
      data = pts,
      mapping = ggplot2::aes(x = .data[["mid"]], y = .data[["med"]] * scale),
      colour = y2col
    ) +
    ggplot2::scale_y_continuous(
      "number of fires",
      sec.axis = ggplot2::sec_axis(
        ~ . / scale,
        name = paste0("median log[fire size] (", size_unit, ")")
      )
    ) +
    ggplot2::labs(x = paste0("log[fire size] (", size_unit, ")"), title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y.right = ggplot2::element_text(colour = y2col),
      axis.text.y.right = ggplot2::element_text(colour = y2col)
    )
}
