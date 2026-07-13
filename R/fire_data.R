## Loaders for the national fire records (NBAC + NFDB polygons, and NFDB points),
## clipped to a study area and harmonised to a common `YEAR` + `SIZE_HA` schema.
## These feed the observed (historical) side of fire-regime summaries.
## The study area may be a path, `sf`, `SpatVector`, or `SpatRaster`, and the year/burned-area column
## names are detected tolerantly across NBAC/NFDB vintages and across projects.

## Pick the first candidate column present, or NA when none match (callers decide
## whether a missing column is fatal). NFDB/NBAC schemas differ between vintages
## and projects (e.g. YEAR vs FIRE_YEAR; ADJ_HA vs POLY_HA vs HECTARES).
.first_col <- function(x, candidates) {
  intersect(candidates, names(x))[1L]
}

## Resolve a study-area argument to a terra object usable for CRS + extent cropping.
## Accepts a path (vector or raster), `sf`/`sfc`, `SpatVector`, or `SpatRaster`.
## A `SpatRaster` is kept as-is (terra::crop() clips a vector to its extent;
## terra::crs() gives its CRS).
.as_study_area <- function(x) {
  if (inherits(x, c("SpatVector", "SpatRaster"))) {
    return(x)
  }
  if (inherits(x, c("sf", "sfc"))) {
    return(terra::vect(x))
  }
  if (is.character(x) && length(x) == 1L) {
    return(tryCatch(terra::vect(x), error = function(e) terra::rast(x)))
  }
  stop("`study_area` must be a file path, sf, SpatVector, or SpatRaster.", call. = FALSE)
}

## Shared body: read shapefile(s), optionally repair invalid geometries, harmonise
## YEAR + SIZE_HA (tolerant columns), filter to fire years + >= `min_size_ha`,
## project + crop to the study area. `size_required = TRUE` (NBAC) errors on a
## missing area column; `FALSE` (NFDB) tolerates it (NFDB always ships SIZE_HA, but
## stay lenient). `repair = TRUE` for polygons (which may be topologically invalid);
## points are always valid, so point loaders pass `repair = FALSE` (skips a needless
## makeValid). Records with a missing (`NA`) size pass the size filter regardless.
.load_fire_vect <- function(
  shp,
  study_area,
  fire_years,
  year_cols,
  size_cols,
  size_required,
  min_size_ha = 1,
  repair = TRUE
) {
  sa <- .as_study_area(study_area)
  p <- lapply(shp, function(x) {
    pp <- withCallingHandlers(terra::vect(x), warning = function(w) {
      if (grepl("Z coordinates ignored", conditionMessage(w))) invokeRestart("muffleWarning")
    })
    if (isTRUE(repair)) {
      pp <- spatialutils::repair_geoms(pp) ## repair invalid geometries (only the invalid subset)
    }
    pp
  }) |>
    tidyterra::bind_spat_rows() ## robust to column differences between multi-year partitions

  year_col <- .first_col(p, year_cols)
  size_col <- .first_col(p, size_cols)
  if (is.na(year_col) || (isTRUE(size_required) && is.na(size_col))) {
    stop(
      sprintf(
        "fire-perimeter shapefile is missing expected year/size columns. Found: %s",
        paste(names(p), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  p <- tidyterra::mutate(
    p,
    YEAR = as.integer(.data[[year_col]]),
    SIZE_HA = if (is.na(size_col)) NA_real_ else as.numeric(.data[[size_col]])
  )
  p <- tidyterra::filter(
    p,
    .data$YEAR %in% !!fire_years,
    is.na(.data$SIZE_HA) | .data$SIZE_HA >= !!min_size_ha
  )
  p <- terra::project(p, terra::crs(sa))
  terra::crop(p, sa) ## clip to the study-area extent
}

#' Load NBAC fire perimeters, harmonised + clipped to a study area
#'
#' Loads National Burned Area Composite (NBAC) fire perimeters and harmonises them
#' to a common `YEAR` + `SIZE_HA` schema: `YEAR` from the NBAC year field and
#' `SIZE_HA` from the adjusted burned area (NBAC's canonical burned-area figure,
#' excluding unburned islands/interior water). Both columns are detected tolerantly
#' (year: `YEAR` or `FIRE_YEAR`; area: `ADJ_HA`, `POLY_HA`, or `HECTARES`) across
#' NBAC vintages. Invalid geometries are repaired with
#' [spatialutils::repair_geoms()] (only the invalid subset is passed to
#' `terra::makeValid()`, for speed), records are filtered to `fire_years` and
#' `SIZE_HA >= min_size_ha`, then projected + cropped to `study_area`.
#'
#' NBAC perimeters are satellite-derived (best-available delineation) and span
#' 1972-present; they are preferred over the NFDB polygon record, whose older
#' perimeters are aerial sketches that overestimate burned area. Use
#' [load_nfdb_polys()] only to backfill years NBAC does not cover.
#'
#' @param nbac_shp Path(s) to the NBAC polygon shapefile(s).
#' @param study_area Study area defining the output CRS + crop extent: a file path
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster` (e.g. a simulation
#'   `flammableMap`).
#' @param fire_years Integer vector of fire years to keep.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept.
#'
#' @returns A `SpatVector` of NBAC perimeters cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @family fire-record loaders
#' @export
load_nbac_polys <- function(nbac_shp, study_area, fire_years, min_size_ha = 1) {
  .load_fire_vect(
    nbac_shp,
    study_area,
    fire_years,
    year_cols = c("YEAR", "FIRE_YEAR"),
    size_cols = c("ADJ_HA", "POLY_HA", "HECTARES"),
    size_required = TRUE,
    min_size_ha = min_size_ha
  )
}

#' Load NFDB fire polygons, harmonised + clipped to a study area
#'
#' Loads National Fire DataBase (NFDB) fire polygons and harmonises them to the
#' same `YEAR` + `SIZE_HA` schema as [load_nbac_polys()]. The year column is
#' detected tolerantly (`YEAR` or `FIRE_YEAR`) and the area column from `SIZE_HA`,
#' `POLY_HA`, or `HECTARES`. Invalid geometries are repaired with
#' [spatialutils::repair_geoms()] (only the invalid subset is passed to
#' `terra::makeValid()`, for speed), records are filtered to `fire_years` and
#' `SIZE_HA >= min_size_ha`, then projected + cropped to `study_area`.
#' The NFDB poly record ships multiple multi-year partitions with differing
#' columns, so pass all their paths together.
#'
#' Prefer NBAC ([load_nbac_polys()]); use NFDB polygons only to backfill years NBAC
#' does not cover.
#'
#' @param nfdb_shp Character vector of NFDB polygon shapefile path(s).
#' @param study_area Study area defining the output CRS + crop extent: a file path
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster`.
#' @param fire_years Integer vector of fire years to keep.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept.
#'
#' @returns A `SpatVector` of NFDB polygons cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @family fire-record loaders
#' @export
load_nfdb_polys <- function(nfdb_shp, study_area, fire_years, min_size_ha = 1) {
  .load_fire_vect(
    nfdb_shp,
    study_area,
    fire_years,
    year_cols = c("YEAR", "FIRE_YEAR"),
    size_cols = c("SIZE_HA", "POLY_HA", "HECTARES"),
    size_required = FALSE,
    min_size_ha = min_size_ha
  )
}

#' Load NFDB fire points, harmonised + clipped to a study area
#'
#' Loads National Fire DataBase (NFDB) fire-point records (the point shapefile of
#' fire locations, e.g. `NFDB_point`) and harmonises them to the same `YEAR` +
#' `SIZE_HA` schema as the perimeter loaders. The year column is detected tolerantly
#' (`YEAR` or `FIRE_YEAR`) and the area column from `SIZE_HA`, `POLY_HA`, or
#' `HECTARES`; records are filtered to `fire_years` and `SIZE_HA >= min_size_ha`,
#' then projected + cropped to `study_area`. Point geometries are always
#' topologically valid, so (unlike the polygon loaders) no geometry repair is
#' performed.
#'
#' NFDB points are fire *locations* rather than mapped burned-area polygons: they
#' carry reported sizes but no perimeter geometry, while covering fires that were
#' never delineated. Use this loader when point locations + reported sizes suffice
#' (e.g. ignition or fire-size-distribution summaries); use the perimeter loaders
#' ([load_nbac_polys()], [load_nfdb_polys()]) when mapped burned-area geometry is
#' required. To retain the small fires that are the point record's main advantage,
#' set `min_size_ha = 0` (the default `1` ha floor matches the perimeter loaders).
#'
#' @param nfdb_shp Character vector of NFDB point shapefile path(s) (typically the
#'   single `NFDB_point` shapefile).
#' @param study_area Study area defining the output CRS + crop extent: a file path
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster`.
#' @param fire_years Integer vector of fire years to keep.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept. Set to `0` to keep all reported fire sizes.
#'
#' @returns A `SpatVector` of NFDB fire points cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @family fire-record loaders
#' @export
load_nfdb_points <- function(nfdb_shp, study_area, fire_years, min_size_ha = 1) {
  .load_fire_vect(
    nfdb_shp,
    study_area,
    fire_years,
    year_cols = c("YEAR", "FIRE_YEAR"),
    size_cols = c("SIZE_HA", "POLY_HA", "HECTARES"),
    size_required = FALSE,
    min_size_ha = min_size_ha,
    repair = FALSE
  )
}
