## Across-replicate five-number summary of `value_col`, pushed down to Arrow compute.
## Arrow's median/quantile are approximate (documented); those expected notes are muffled.
.dist_summary <- function(ds, value_col, id_cols, count_name) {
  vsym <- rlang::sym(value_col)
  out <- withCallingHandlers(
    ds |>
      dplyr::group_by(!!!rlang::syms(id_cols)) |>
      dplyr::summarise(
        .n = dplyr::n(),
        mean = mean(!!vsym, na.rm = TRUE),
        sd = stats::sd(!!vsym, na.rm = TRUE),
        min = min(!!vsym, na.rm = TRUE),
        q25 = stats::quantile(!!vsym, 0.25, na.rm = TRUE),
        median = stats::median(!!vsym, na.rm = TRUE),
        q75 = stats::quantile(!!vsym, 0.75, na.rm = TRUE),
        max = max(!!vsym, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::collect(),
    warning = function(w) {
      if (grepl("approximate (median|quantile)", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  out <- as.data.frame(out)
  names(out)[names(out) == ".n"] <- count_name
  out$se <- out$sd / sqrt(out[[count_name]])
  out$ci <- out$se * stats::qt(0.975, pmax(out[[count_name]] - 1L, 1L))
  out
}

.resolve_dataset <- function(x, opener) {
  if (inherits(x, c("Dataset", "arrow_dplyr_query"))) x else opener(x)
}

#' Summarize the across-replicate fire-size distribution
#'
#' Opens the per-replicate fire-event parquet(s) lazily and pushes the reduction
#' down to Arrow compute, so the per-event rows are never all held in memory at
#' once. All fire events in each identifier group (pooled across replicates and,
#' by default, across years) contribute to the distribution: the count and the
#' five-number summary (`min`, `q25`, `median`, `q75`, `max`) of `size_col` plus
#' the mean, sd, and standard error / 95% confidence half-width. The quantiles and
#' median are Arrow-approximate; `min`/`max`/`mean`/`sd`/`n_fires` are exact.
#'
#' @param x Parquet paths / roots (see [open_burn_dataset()]) or an Arrow
#'   `Dataset` / query.
#' @param size_col Name of the fire-size column (default `"size"`).
#' @param id_cols Grouping columns; `NULL` (default) auto-detects those present
#'   among `simArea`, `studyArea`, `scenario`, `region`.
#'
#' @return A `data.frame` with the id columns plus `n_fires`, `mean`, `sd`, `min`,
#'   `q25`, `median`, `q75`, `max`, `se`, `ci`; zero rows if there is no data.
#'
#' @export
summarize_fire_sizes <- function(x, size_col = "size", id_cols = NULL) {
  ds <- .resolve_dataset(x, open_burn_dataset)
  if (is.null(ds)) {
    return(data.frame())
  }
  cols <- names(ds)
  if (!size_col %in% cols) {
    stop(
      "summarize_fire_sizes(): size column '",
      size_col,
      "' not in dataset (have: ",
      paste(cols, collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  if (is.null(id_cols)) {
    id_cols <- intersect(c("simArea", "studyArea", "scenario", "region"), cols)
  }
  id_cols <- setdiff(id_cols, size_col)
  .dist_summary(ds, size_col, id_cols, "n_fires")
}

#' Summarize across-replicate annual area burned
#'
#' Like [summarize_fire_sizes()] but for a per-replicate annual-area-burned table:
#' for each identifier group (by default one per time step and region) it returns
#' the replicate count and the across-replicate distribution of `area_col`.
#'
#' @param x Parquet paths / roots (see [open_burn_dataset()]) or an Arrow
#'   `Dataset` / query.
#' @param area_col Name of the area-burned column (default `"area_ha"`).
#' @param id_cols Grouping columns; `NULL` (default) auto-detects those present
#'   among `simArea`, `studyArea`, `scenario`, `time`, `year`, `region`.
#'
#' @return A `data.frame` with the id columns plus `n_reps`, `mean`, `sd`, `min`,
#'   `q25`, `median`, `q75`, `max`, `se`, `ci`; zero rows if there is no data.
#'
#' @export
summarize_annual_area <- function(x, area_col = "area_ha", id_cols = NULL) {
  ds <- .resolve_dataset(x, open_burn_dataset)
  if (is.null(ds)) {
    return(data.frame())
  }
  cols <- names(ds)
  if (!area_col %in% cols) {
    stop(
      "summarize_annual_area(): area column '",
      area_col,
      "' not in dataset (have: ",
      paste(cols, collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  if (is.null(id_cols)) {
    id_cols <- intersect(c("simArea", "studyArea", "scenario", "time", "year", "region"), cols)
  }
  id_cols <- setdiff(id_cols, area_col)
  .dist_summary(ds, area_col, id_cols, "n_reps")
}

#' Fire cycle from a per-replicate annual-area-burned table
#'
#' The fire cycle is the number of years for an area equal to the whole landscape
#' to burn: `landscape_ha / mean(area_col)`, where the mean is taken over all rows
#' (pooling replicates and years) within each identifier group. Pre-filter the
#' input (e.g. to a single region) before calling if needed.
#'
#' @param x Parquet paths / roots (see [open_burn_dataset()]) or an Arrow
#'   `Dataset` / query.
#' @param area_col Name of the annual-area-burned column (default `"area_ha"`).
#' @param landscape_ha Total flammable/active landscape area, in hectares.
#' @param id_cols Grouping columns; `NULL` (default) auto-detects those present
#'   among `simArea`, `studyArea`, `scenario`.
#'
#' @return A `data.frame` with the id columns plus `mean_annual_area` and
#'   `fire_cycle_yr`; zero rows if there is no data.
#'
#' @export
fire_cycle <- function(x, area_col = "area_ha", landscape_ha, id_cols = NULL) {
  ds <- .resolve_dataset(x, open_burn_dataset)
  if (is.null(ds)) {
    return(data.frame())
  }
  cols <- names(ds)
  if (!area_col %in% cols) {
    stop(
      "fire_cycle(): area column '",
      area_col,
      "' not in dataset (have: ",
      paste(cols, collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  if (is.null(id_cols)) {
    id_cols <- intersect(c("simArea", "studyArea", "scenario"), cols)
  }
  id_cols <- setdiff(id_cols, area_col)
  asym <- rlang::sym(area_col)
  out <- ds |>
    dplyr::group_by(!!!rlang::syms(id_cols)) |>
    dplyr::summarise(mean_annual_area = mean(!!asym, na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  out <- as.data.frame(out)
  out$fire_cycle_yr <- landscape_ha / out$mean_annual_area
  out
}
