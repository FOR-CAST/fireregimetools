## Fire causes and prescribed burns in the national fire records.
##
## The value tables below follow the records' own documentation and the values found in them:
## - NBAC polygons: `FIRECAUS` is "Natural", "Human" or "Undetermined"; `PRESCRIBED` is "true" or blank
##   (NBAC 1972-2025, release 2026-05-13).
## - NFDB polygons and points: `CAUSE` is N (natural, e.g. lightning; L in older releases), H (human),
##   H-PB (prescribed burn), U (unknown), Re (reburn: the original cause was not recorded), or "n/a"
##   (NFDB polygons 2021-07-07). From the 2026 point release, `CAUSE` holds only N, H and U, and the
##   H-PB and RE codes are in `CAUSE2`. `PRESCRIBED` is "PB", "1", "No", "0" or blank on points, and
##   `FIRE_TYPE` is "PB" for some prescribed burns.

## lower-case cause value -> cause group, by record; a blank or missing value is "unknown"
.fire_cause_values <- list(
  nbac = c(natural = "natural", human = "human", undetermined = "unknown", unknown = "unknown"),
  nfdb = c(
    n = "natural",
    l = "natural",
    h = "human",
    `h-pb` = "human",
    u = "unknown",
    `n/a` = "unknown",
    re = "other"
  )
)

.fire_cause_groups <- c("natural", "unknown", "human", "other")

## the default cause column of each record
.fire_cause_col <- function(source) if (source == "nbac") "FIRECAUS" else "CAUSE"

## `x` as lower-case, trimmed text, with NA as ""
.fire_value <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

#' Cause group of each fire record
#'
#' Classifies the cause each NBAC or NFDB record gives into one of four groups: `"natural"`
#' (NBAC `Natural`; NFDB `N`, or `L` in older releases), `"human"` (NBAC `Human`; NFDB `H` and `H-PB`,
#' a prescribed burn being human-caused), `"unknown"` (NBAC `Undetermined`; NFDB `U` and `n/a`; and a
#' blank or missing cause in either), and `"other"` (NFDB `Re`, a reburn whose original cause was not
#' recorded). Values are matched ignoring case and surrounding spaces.
#'
#' @param records `SpatVector`, `sf` or `data.frame` of fire records.
#' @param source `"nbac"` or `"nfdb"`: which record `records` come from.
#' @param cause_col Column holding the cause. Default `"FIRECAUS"` for NBAC, `"CAUSE"` for NFDB.
#'
#' @return A character vector, one cause group per record.
#'   A value not listed above is an error, rather than a guess.
#'
#' @seealso [filter_fire_records()], [is_prescribed_burn()]
#' @export
fire_cause_group <- function(records, source = c("nbac", "nfdb"), cause_col = NULL) {
  source <- match.arg(source)
  cause_col <- cause_col %||% .fire_cause_col(source)
  d <- as.data.frame(records)
  if (!cause_col %in% names(d)) {
    stop(
      sprintf("the fire records have no `%s` column, so their causes are not known", cause_col),
      call. = FALSE
    )
  }
  values <- .fire_cause_values[[source]]
  v <- .fire_value(d[[cause_col]])
  group <- unname(values[match(v, names(values))])
  group[v == ""] <- "unknown"
  if (anyNA(group)) {
    stop(
      sprintf(
        "unrecognised %s `%s` value(s): %s",
        toupper(source),
        cause_col,
        paste(unique(d[[cause_col]][is.na(group)]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  group
}

#' Which fire records are prescribed burns
#'
#' A record is a prescribed burn when any of these holds (ignoring case): its `PRESCRIBED` is
#' `"true"`, `"t"`, `"yes"`, `"y"`, `"1"` or `"PB"`; its `FIRE_TYPE` is `"PB"` or `"Prescribed Burn"`;
#' or its `CAUSE` or `CAUSE2` is `"H-PB"`. A column the records lack is not consulted, and a blank or
#' missing value, `"No"`, `"0"` or `"false"` does not mark a prescribed burn.
#'
#' @inheritParams fire_cause_group
#'
#' @return A logical vector, one per record: `TRUE` for a prescribed burn.
#'
#' @seealso [filter_fire_records()]
#' @export
is_prescribed_burn <- function(records) {
  d <- as.data.frame(records)
  has <- function(col, values) {
    if (col %in% names(d)) .fire_value(d[[col]]) %in% values else rep(FALSE, nrow(d))
  }
  has("PRESCRIBED", c("true", "t", "yes", "y", "1", "pb")) |
    has("FIRE_TYPE", c("pb", "prescribed burn")) |
    has("CAUSE", "h-pb") |
    has("CAUSE2", "h-pb")
}

#' Keep the fire records of some causes, without prescribed burns
#'
#' Selects the fire records a fire model is compared with: by default every record that is not a
#' prescribed burn ([is_prescribed_burn()]), and with `causes`, only those whose cause group
#' ([fire_cause_group()]) is one of them.
#'
#' @inheritParams fire_cause_group
#' @param causes `NULL` (default) keeps every cause. Otherwise cause groups to keep, any of
#'   `"natural"`, `"unknown"`, `"human"` and `"other"`.
#' @param drop_prescribed Drop prescribed burns (default `TRUE`), whatever `causes` says.
#'
#' @return `records` without the rows left out.
#'
#' @seealso [load_nbac_polys()], [load_nfdb_polys()] and [load_nfdb_points()], which keep the cause
#'   and prescribed-burn columns of the records they read.
#' @export
filter_fire_records <- function(
  records,
  causes = NULL,
  source = c("nbac", "nfdb"),
  cause_col = NULL,
  drop_prescribed = TRUE
) {
  source <- match.arg(source)
  bad <- setdiff(causes, .fire_cause_groups)
  if (length(bad)) {
    stop(
      sprintf(
        "unknown fire cause group(s): %s (use %s)",
        paste(bad, collapse = ", "),
        paste(.fire_cause_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  keep <- rep(TRUE, nrow(records))
  if (isTRUE(drop_prescribed)) {
    keep <- !is_prescribed_burn(records)
  }
  if (!is.null(causes)) {
    keep <- keep & fire_cause_group(records, source, cause_col) %in% causes
  }
  records[which(keep), ]
}
