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

## Muffle terra's "Z coordinates ignored" warning: several NBAC/NFDB vintages carry Z values and
## trigger it on every read, and it says nothing a caller can act on.
.muffle_z <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Z coordinates ignored", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

## Records with no CRS cannot be projected + cropped, and terra::project() rejects them with a
## message that names neither the argument nor the file at fault.
.stop_if_no_crs <- function(crs_wkt, files) {
  if (!nzchar(crs_wkt)) {
    stop(
      "the fire records have no CRS (is a .prj missing?): ",
      paste(basename(unlist(files)), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## A conservative superset of `sa`'s footprint, expressed in `src_crs`, for pushing a spatial filter
## down to the GDAL/OGR read. Returns NULL when no trustworthy extent can be derived, which makes
## the caller read the whole source -- slow, but never wrong.
##
## Why a PADDED BOUNDING BOX, rather than the study-area geometry:
##
##  * the authoritative selection is still made downstream by terra::crop(), in the study area's own
##    CRS. This filter's only job is to keep the read small, so it must err towards selecting too
##    much: the surplus is discarded by crop() anyway, while a record missed here is missed for
##    good. terra::project() moves a polygon's VERTICES only, so a study-area outline projected into
##    the source CRS has straight chords where the true boundary curves -- a filter built from it can
##    fall INSIDE the real footprint and drop an edge record that crop() would have kept.
##  * a `SpatRaster` study area selects by extent by documented contract, so a geometry-exact
##    pre-filter would quietly change that too.
##  * terra's pushdown is bbox-based regardless: a `filter=` polygon whose footprint is L-shaped
##    still returns the records sitting in its notch.
##
## How the box is built: the study area's extent is widened, its outline densified along STRAIGHT
## lines in the study area's CRS (`flat = TRUE`) -- which is how crop() draws its edges -- and the
## box taken around that outline's projected vertices. The corners alone are not enough, because a
## reprojected edge curves: a lon/lat study area's parallels bow in BC Albers, so a record just inside
## the southern edge falls outside the extent of the projected corners. Nor is terra's own
## project(<SpatExtent>), which densifies along great circles; those leave a parallel edge, and missed
## exactly such a record (#1, which contributed this construction and its regression test).
##
## An outline's projected box bounds the whole study area wherever the projection is continuous,
## which leaves two ways it can fail, and both are detected rather than padded over:
##
##  * project() silently DROPS points it cannot convert (leaving an empty geometry rather than a
##    non-finite coordinate), so an outline that does not come back whole means the study area
##    reaches outside the source CRS's domain -- and a box around only the survivors is too small.
##  * a lon/lat source is discontinuous at the poles and the antimeridian. A study area containing a
##    pole maps its outline to a ring of latitudes that stops short of it, and one straddling the
##    antimeridian maps it to both ends of the longitude range. Either way the outline's longitudes
##    wrap (span > 180 degrees), and no box around them can be trusted.
##
## `pad_frac` slack is added twice, once in each CRS, so the box covers a *buffered* study area. That
## guards the one vertex-only-projection effect the outline cannot: a source geometry lying outside
## the study area can, once its own vertices are projected and re-joined by straight chords, cut a
## corner into it. Only a source segment longer than the pad could do so, which no perimeter in these
## records comes close to. `pad_frac = 0` is allowed, so the outline can be tested on its own.
.prefilter_extent <- function(sa, src_crs, pad_frac = 0.05, interval_frac = 0.01) {
  if (!nzchar(src_crs)) {
    return(NULL)
  }
  e <- terra::ext(sa)
  ev <- as.vector(e)
  if (!all(is.finite(ev)) || ev[["xmax"]] <= ev[["xmin"]] || ev[["ymax"]] <= ev[["ymin"]]) {
    return(NULL)
  }
  span <- max(ev[["xmax"]] - ev[["xmin"]], ev[["ymax"]] - ev[["ymin"]])
  ## Do NOT change this to `flat = FALSE`. Beyond densifying along great circles (the wrong edge, see
  ## above), `flat = FALSE` changes the UNITS of `interval`: `flat = TRUE` reads it in the study
  ## area's CRS units, but a geodesic densify reads it in METRES. On a lon/lat study area this
  ## degree-sized interval then means centimetres -- an 8 x 6 degree study area densified at 0.08 gets
  ## 351 vertices with `flat = TRUE` and 30,094,846 with `flat = FALSE`, and projecting those peaks at
  ## 27.6 GB (measured) -- enough to OOM a workstation. A geodesic outline would need its interval
  ## converted to metres, and would still be the wrong edge.
  outline <- terra::as.polygons(terra::extend(e, pad_frac * span), crs = terra::crs(sa)) |>
    terra::densify(interval_frac * span, flat = TRUE) |>
    terra::as.points()
  ## this projects an internal scaffold, not user data, so its warnings (out-of-domain points, datum
  ## shifts) are noise; a projection that actually fails is caught by the NULL returns instead
  xy <- tryCatch(
    suppressWarnings(terra::crds(terra::project(outline, src_crs))),
    error = function(e) NULL
  )
  if (is.null(xy) || nrow(xy) != nrow(terra::crds(outline)) || !all(is.finite(xy))) {
    return(NULL)
  }
  xr <- range(xy[, 1L])
  yr <- range(xy[, 2L])
  if (isTRUE(terra::is.lonlat(src_crs)) && xr[2L] - xr[1L] > 180) {
    return(NULL)
  }
  span_src <- max(xr[2L] - xr[1L], yr[2L] - yr[1L])
  if (!is.finite(span_src) || span_src <= 0) {
    return(NULL)
  }
  terra::extend(terra::ext(xr[1L], xr[2L], yr[1L], yr[2L]), pad_frac * span_src)
}

## Read one fire-record source, pushing the study-area spatial filter down to the GDAL/OGR read so
## that only candidate features are materialised in R. The national records are far larger than any
## study area (the NBAC shapefile alone is ~1.9 GB for the whole of Canada), and everything read here
## is then repaired, projected and cropped -- so filtering at the read is what keeps that work
## proportional to the study area rather than to the country.
##
## `terra::vect(proxy = TRUE)` opens the source without reading geometry, which also lets a missing
## source CRS be reported before gigabytes have been read. Every step falls back to reading the whole
## source, so the pushdown can only shrink the read, never change the result.
##
## Only the SPATIAL filter is pushed down, not `fire_years`/`min_size_ha`: the geometry is what costs
## the memory, the attribute filters are cheap once the spatial filter has done its work, and an
## OGR-SQL `where` clause would have to second-guess each vintage's column types (the loaders already
## tolerate a `YEAR` that needs coercing) for no measurable gain.
##
## spatialutils::read_vector_aoi() does the same pushdown, but it refines the read with
## terra::relate() against the AOI, which makes the SOURCE-CRS geometry the authoritative selection.
## Here it must not be: the study area's own CRS is where records are selected (see
## .prefilter_extent()), and a `SpatRaster` study area selects by extent by documented contract.
.read_fire_source <- function(x, sa, prefilter = TRUE) {
  ## a fallback is correct but can turn seconds into many minutes on the national records, so say so
  ## rather than leave the caller wondering where the time went
  read_all <- function(why = NULL) {
    if (!is.null(why)) {
      message("no spatial pre-filter for ", basename(x), " (", why, "); reading the whole source")
    }
    .muffle_z(terra::vect(x))
  }
  if (!isTRUE(prefilter)) {
    return(read_all())
  }
  px <- tryCatch(terra::vect(x, proxy = TRUE), error = function(e) NULL)
  if (is.null(px)) {
    return(read_all("it cannot be opened for a filtered read"))
  }
  .stop_if_no_crs(terra::crs(px), x)
  e <- .prefilter_extent(sa, terra::crs(px))
  if (is.null(e)) {
    return(read_all("the study area has no usable footprint in its CRS"))
  }
  out <- tryCatch(.muffle_z(terra::query(px, extent = e)), error = function(err) NULL)
  if (is.null(out)) read_all("the filtered read failed") else out
}

## Bind the per-source reads into one SpatVector, keeping the attribute schema when nothing was
## selected. tidyterra::bind_spat_rows() is what makes differing columns across the NFDB polygon
## record's multi-year partitions bindable, but handed nothing but empty parts it returns a
## single placeholder `first_empty` column instead of the schema -- and a study area that selects no
## records is ordinary (a pushed-down read returns none, and so does a study area outside the
## records' footprint), so callers would then be told their file was missing its year/size columns.
## Empty parts carry nothing to bind, so drop them; if that leaves nothing, hand back the first read
## as the (zero-row) schema carrier.
.bind_fire_parts <- function(parts) {
  keep <- parts[vapply(parts, function(x) nrow(x) > 0L, logical(1L))]
  if (!length(keep)) {
    return(parts[[1L]])
  }
  tidyterra::bind_spat_rows(keep)
}

## Shared body: read shapefile(s), optionally repair invalid geometries, harmonise
## YEAR + SIZE_HA (tolerant columns), filter to fire years + >= `min_size_ha`,
## project + crop to the study area. `size_required = TRUE` (NBAC) errors on a
## missing area column; `FALSE` (NFDB) tolerates it (NFDB always ships SIZE_HA, but
## stay lenient). `repair = TRUE` for polygons (which may be topologically invalid);
## points are always valid, so point loaders pass `repair = FALSE` (skips a needless
## makeValid). Records with a missing (`NA`) size pass the size filter regardless.
## `prefilter = FALSE` reads each source whole instead of pushing the study-area filter down to
## GDAL: the same records, read the slow way -- an escape hatch, and the reference path the tests
## compare the pushdown against.
.load_fire_vect <- function(
  shp,
  study_area,
  fire_years = NULL,
  year_cols,
  size_cols,
  size_required,
  min_size_ha = 1,
  repair = TRUE,
  prefilter = getOption("fireregimetools.prefilter", TRUE)
) {
  sa <- .as_study_area(study_area)
  ## Check the study area's CRS up front. terra::project() does reject an empty target CRS, but only
  ## after every shapefile has been read and repaired -- minutes and ~1 GB of I/O later -- and it
  ## says "[project] output crs is not valid" without naming the argument at fault.
  if (!nzchar(terra::crs(sa))) {
    stop(
      "`study_area` has no CRS; one is required to project + crop the fire records.",
      call. = FALSE
    )
  }
  parts <- lapply(shp, function(x) {
    ## the study-area filter is pushed down to the read, so the repair below (and the projection and
    ## crop further down) see only the records near the study area rather than the whole country.
    ## Selecting before repairing rather than after is safe: the filter tests each record's STORED
    ## bounding box, and makeValid() derives its output from the input's own linework, so a repaired
    ## geometry never reaches outside the box its broken form was selected by.
    pp <- .read_fire_source(x, sa, prefilter = prefilter)
    if (isTRUE(repair)) {
      pp <- spatialutils::repair_geoms(pp) ## repair invalid geometries (only the invalid subset)
    }
    pp
  })
  p <- .bind_fire_parts(parts)

  ## likewise for records with no CRS: terra::project() rejects those too, with a message that says
  ## nothing about which file is at fault. The pushdown path reports this from the proxy, before any
  ## geometry is read; this backstops the whole-source reads, which do not consult one.
  .stop_if_no_crs(terra::crs(p), shp)

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
  if (!is.null(fire_years)) {
    p <- tidyterra::filter(p, .data$YEAR %in% !!fire_years)
  }
  p <- tidyterra::filter(p, is.na(.data$SIZE_HA) | .data$SIZE_HA >= !!min_size_ha)
  p <- terra::project(p, terra::crs(sa))
  out <- terra::crop(p, sa) ## geometry for a vector study area; extent for a SpatRaster
  ## terra::crop() returns an attribute-LESS SpatVector (ncol 0) when nothing survives, so callers
  ## that reach for `$YEAR` get NULL rather than an empty vector. Hand back an empty subset of the
  ## harmonised input instead, which keeps the schema.
  if (nrow(out) == 0L) p[0, ] else out
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
#' Perimeters straddling the study-area boundary are clipped, so their geometry is
#' truncated while `SIZE_HA` keeps the full reported fire size. Sum `SIZE_HA` for
#' "fires that reached this study area"; use [terra::expanse()] for "area burned
#' inside it" -- the two differ for every edge fire.
#'
#' NBAC perimeters are satellite-derived (best-available delineation) and span
#' 1972-present; they are preferred over the NFDB polygon record, whose older
#' perimeters are aerial sketches that overestimate burned area. Use
#' [load_nfdb_polys()] only to backfill years NBAC does not cover.
#'
#' @param nbac_shp Path(s) to the NBAC polygon shapefile(s).
#' @param study_area Study area defining the output CRS + crop extent: a file path
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster` (e.g. a simulation
#'   `flammableMap`). It must carry a CRS. A **vector** study area selects records by
#'   its geometry; a **`SpatRaster`** selects by its extent, and its `NA` cells do not
#'   narrow that selection -- so an irregular study area passed as a raster also
#'   returns records from the bounding box around it.
#' @param fire_years Integer vector of fire years to keep; `NULL` (default) keeps all years.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept.
#'
#' @returns A `SpatVector` of NBAC perimeters cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @section Spatial pre-filtering:
#' The study-area filter is pushed down to the GDAL/OGR read, so only records near
#' the study area are materialised in R: the source is opened as a
#' [terra::vect()] *proxy* (no geometry read), the study area's footprint is
#' projected into the records' own CRS as a padded bounding box, and only the
#' features inside it are read. Geometry repair, projection and the final crop
#' then see a study-area-sized subset rather than the whole national record, which
#' is what keeps the loaders' time and memory proportional to the study area.
#'
#' The box is deliberately generous -- the exact selection is still made by
#' `terra::crop()` in the study area's CRS, and every step of the pushdown falls
#' back to reading the source whole -- so results do not depend on it. Set
#' `options(fireregimetools.prefilter = FALSE)` to read the sources whole
#' regardless; the records returned are the same, they just take longer to arrive.
#'
#' @family fire-record loaders
#' @export
load_nbac_polys <- function(nbac_shp, study_area, fire_years = NULL, min_size_ha = 1) {
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
#' Perimeters straddling the study-area boundary are clipped, so their geometry is
#' truncated while `SIZE_HA` keeps the full reported fire size. Sum `SIZE_HA` for
#' "fires that reached this study area"; use [terra::expanse()] for "area burned
#' inside it" -- the two differ for every edge fire.
#'
#' Prefer NBAC ([load_nbac_polys()]); use NFDB polygons only to backfill years NBAC
#' does not cover.
#'
#' @param nfdb_shp Character vector of NFDB polygon shapefile path(s).
#' @param study_area Study area defining the output CRS + crop extent: a file path
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster` (e.g. a simulation
#'   `flammableMap`). It must carry a CRS. A **vector** study area selects records by
#'   its geometry; a **`SpatRaster`** selects by its extent, and its `NA` cells do not
#'   narrow that selection -- so an irregular study area passed as a raster also
#'   returns records from the bounding box around it.
#' @param fire_years Integer vector of fire years to keep; `NULL` (default) keeps all years.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept.
#'
#' @returns A `SpatVector` of NFDB polygons cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @family fire-record loaders
#' @inheritSection load_nbac_polys Spatial pre-filtering
#' @export
load_nfdb_polys <- function(nfdb_shp, study_area, fire_years = NULL, min_size_ha = 1) {
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
#'   (vector or raster), `sf`, `SpatVector`, or `SpatRaster` (e.g. a simulation
#'   `flammableMap`). It must carry a CRS. A **vector** study area selects records by
#'   its geometry; a **`SpatRaster`** selects by its extent, and its `NA` cells do not
#'   narrow that selection -- so an irregular study area passed as a raster also
#'   returns records from the bounding box around it.
#' @param fire_years Integer vector of fire years to keep; `NULL` (default) keeps all years.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records
#'   with a smaller reported `SIZE_HA` are dropped, while records with a missing
#'   (`NA`) size are always kept. Set to `0` to keep all reported fire sizes.
#'
#' @returns A `SpatVector` of NFDB fire points cropped to `study_area`, in its CRS,
#'   carrying harmonised integer `YEAR` + numeric `SIZE_HA` columns.
#'
#' @family fire-record loaders
#' @inheritSection load_nbac_polys Spatial pre-filtering
#' @export
load_nfdb_points <- function(nfdb_shp, study_area, fire_years = NULL, min_size_ha = 1) {
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

## --- Archive acquisition -------------------------------------------------------------------------
##
## The national fire archives are large (NFDB points ~42 MB, NFDB polygons ~780 MB, NBAC ~1.2 GB) and
## are served as plain HTTP downloads, which fail in two SILENT ways this code defends against:
##
##   1. a transfer that exceeds R's default 60 s `download.file` timeout is truncated, not errored; and
##   2. an interrupted extraction leaves a SHORT file behind, which an existence-only cache check then
##      accepts on every later run, permanently.
##
## Both bit us together in 2026: the 1.88 GB NBAC `.shp` sat at 567 MB, GDAL logged 74,178 read
## errors, and the reader still returned the full feature count (the `.shx` index was intact), so a
## whole set of historic fire summaries was built from ~30% of the record with nothing failing.
##
## So: the timeout is raised, the download is staged through a `.part` file and renamed only on
## success, the extraction is staged in a private directory and verified against the archive manifest
## before being published into `dest`, and a stamp file records that a verified extraction happened.

## Marker written only after an extraction has been verified complete.
.fire_archive_stamp <- function(dest, url) {
  file.path(dest, paste0(".", basename(url), ".complete"))
}

## Extracted files matching `pattern`, sorted. `list.files()` skips dotfiles, so the stamp, lock and
## staging directories used below are invisible here by construction. The top level of `dest` is
## preferred and subdirectories are searched only if nothing matches there: `dest` is often a shared
## inputs directory, and a recursive match can otherwise reach an older copy of the same record
## nested somewhere under it -- while an archive that extracts into a folder of its own still works.
.fire_archive_files <- function(dest, pattern) {
  found <- list.files(dest, pattern = pattern, full.names = TRUE)
  if (!length(found)) {
    found <- list.files(dest, pattern = pattern, full.names = TRUE, recursive = TRUE)
  }
  sort(found)
}

## Read a zip's central directory (member `Name` + `Length`). The central directory lives at the END
## of the file, so a truncated archive cannot be listed -- making this an integrity check on the zip
## itself, not just a manifest. `utils::unzip(list = TRUE)` is the primary reader, with libarchive as
## a fallback for zip64 archives it cannot parse; but `archive::archive()` reports every size as 0
## when `options(encoding = "UTF-8")` is set, so an all-zero manifest is treated as unusable rather
## than as "every member is empty". Returns NULL when no usable manifest can be read.
.fire_archive_manifest <- function(zip) {
  m <- tryCatch(suppressWarnings(utils::unzip(zip, list = TRUE)), error = function(e) NULL)
  if (is.null(m) && requireNamespace("archive", quietly = TRUE)) {
    m <- tryCatch(
      {
        a <- archive::archive(zip)
        data.frame(Name = as.character(a$path), Length = as.numeric(a$size))
      },
      error = function(e) NULL
    )
  }
  if (!is.null(m) && (!nrow(m) || all(m$Length == 0))) {
    m <- NULL
  }
  m
}

## Members of `manifest` that are missing from `dir` or shorter than the archive says they should be
## (empty when the extraction is intact, or when there is no manifest to check against).
.fire_archive_shortfall <- function(dir, manifest) {
  if (is.null(manifest)) {
    return(character(0))
  }
  keep <- !grepl("/$", manifest$Name)
  nms <- manifest$Name[keep]
  want <- as.numeric(manifest$Length[keep])
  got <- file.size(file.path(dir, nms))
  nms[is.na(got) | got != want]
}

## Already-extracted files in `dest`, but only when the extraction is trustworthy: this function
## verified it earlier (the stamp), or the archive is still present to verify against now. "A file
## with the right name exists" is exactly the check that accepted a 30%-complete NBAC shapefile, so
## it is not sufficient on its own -- except for files a user placed in `dest` by hand, where there is
## no archive to check against and honouring them is the point; those are accepted with a note.
.fire_archive_cached <- function(dest, url, pattern) {
  found <- .fire_archive_files(dest, pattern)
  if (!length(found)) {
    return(character(0))
  }
  stamp <- .fire_archive_stamp(dest, url)
  if (file.exists(stamp)) {
    return(found)
  }
  zip <- file.path(dest, basename(url))
  if (!file.exists(zip)) {
    message(
      "using existing files in ",
      dest,
      " matching '",
      pattern,
      "' (no archive present to verify them against)"
    )
    return(found)
  }
  manifest <- .fire_archive_manifest(zip)
  if (is.null(manifest) || length(.fire_archive_shortfall(dest, manifest))) {
    return(character(0)) ## incomplete (or unverifiable) -- re-extract below
  }
  try(writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), stamp), silent = TRUE)
  found
}

## Move `from` onto `to`, replacing whatever is there. Rename is atomic within a filesystem, but
## Windows refuses to rename onto an existing file, and a handle left open by a failed read of that
## file (a truncated archive that would not list, say) can outlive the unlink() by a moment. So clear
## the target, retry briefly, and fall back to a copy.
.replace_file <- function(from, to) {
  for (i in seq_len(5L)) {
    unlink(to)
    if (!file.exists(to) && isTRUE(suppressWarnings(file.rename(from, to)))) {
      return(TRUE)
    }
    Sys.sleep(0.2)
  }
  if (isTRUE(suppressWarnings(file.copy(from, to, overwrite = TRUE)))) {
    unlink(from)
    return(TRUE)
  }
  FALSE
}

## Download `url` to `zip` unless a readable (i.e. non-truncated) copy is already there. Staged
## through a process-private `.part` file and renamed only on success, so an interrupted transfer is
## never mistaken for a complete one; the timeout is raised for the multi-hundred-MB archives and R's
## truncation warning is promoted to an error.
.fire_archive_download <- function(url, zip) {
  if (file.exists(zip)) {
    if (!is.null(.fire_archive_manifest(zip))) {
      return(zip)
    }
    message("existing archive is unreadable or incomplete; re-downloading: ", zip)
    unlink(zip)
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(3600L, old_timeout))

  part <- sprintf("%s.part-%d", zip, Sys.getpid())
  ok <- FALSE
  on.exit(if (!ok) unlink(part), add = TRUE)
  status <- withCallingHandlers(
    utils::download.file(url, part, mode = "wb"),
    warning = function(w) {
      if (grepl("downloaded length", conditionMessage(w), fixed = TRUE)) {
        stop("download of ", url, " was truncated: ", conditionMessage(w), call. = FALSE)
      }
    }
  )
  if (!identical(as.integer(status), 0L)) {
    stop("download of ", url, " failed with status ", status, call. = FALSE)
  }
  ## a download that produced nothing (or an empty file) must not be reported as a rename problem
  if (!file.exists(part) || file.size(part) == 0) {
    stop("download of ", url, " produced no data", call. = FALSE)
  }
  if (!.replace_file(part, zip)) {
    stop("unable to move the downloaded archive into place at ", zip, call. = FALSE)
  }
  ok <- TRUE
  zip
}

## Extract with libarchive when available (it handles zip64 and members that R's internal unzip
## cannot), else `utils::unzip()`, whose warnings are promoted to errors so a failed extraction is
## never silently accepted. Paths are absolute because `archive_extract()` changes directory.
.fire_archive_extract <- function(zip, exdir) {
  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
  exdir <- normalizePath(exdir, mustWork = TRUE)
  zip <- normalizePath(zip, mustWork = TRUE)
  if (requireNamespace("archive", quietly = TRUE)) {
    extracted <- tryCatch(
      {
        archive::archive_extract(zip, dir = exdir)
        TRUE
      },
      error = function(e) FALSE
    )
    if (extracted) {
      return(invisible(exdir))
    }
  }
  withCallingHandlers(utils::unzip(zip, exdir = exdir), warning = function(w) {
    stop("extraction of ", basename(zip), " failed: ", conditionMessage(w), call. = FALSE)
  })
  invisible(exdir)
}

## Move a verified extraction from `staging` into `dest`. Rename is atomic within a filesystem
## (`staging` is a subdirectory of `dest`, so it always is), which is what makes publishing safe
## while another process is reading `dest`.
.fire_archive_publish <- function(staging, dest) {
  rel <- list.files(staging, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  for (f in rel) {
    src <- file.path(staging, f)
    tgt <- file.path(dest, f)
    dir.create(dirname(tgt), showWarnings = FALSE, recursive = TRUE)
    if (!.replace_file(src, tgt)) {
      stop("unable to move extracted file into place: ", tgt, call. = FALSE)
    }
  }
  invisible(dest)
}

## Acquire a national fire archive: return ALL extracted files matching `pattern` in `dest` (the NFDB
## polygon record ships several multi-year partitions, which the loaders bind together), downloading
## and extracting `url` only when a verified copy is not already there.
##
## `dest` may be shared by concurrent workers, so one of them takes a lock (an atomically created
## directory) and the others wait rather than each pulling a multi-hundred-MB copy. A worker that
## dies leaves its lock behind, so the wait is bounded by `lock_timeout` seconds; on timeout we fetch
## too, which is safe because extraction is staged privately and published file-by-file.
.download_fire_archive <- function(url, dest, pattern, lock_timeout = 3600) {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  found <- .fire_archive_cached(dest, url, pattern)
  if (length(found)) {
    return(found)
  }

  lock <- file.path(dest, paste0(".", basename(url), ".lock"))
  if (dir.create(lock, showWarnings = FALSE)) {
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  } else {
    message("waiting for another process to fetch ", basename(url), " ...")
    waited <- 0
    while (dir.exists(lock) && waited < lock_timeout) {
      Sys.sleep(5)
      waited <- waited + 5
    }
    found <- .fire_archive_cached(dest, url, pattern)
    if (length(found)) {
      return(found)
    }
  }

  zip <- .fire_archive_download(url, file.path(dest, basename(url)))
  manifest <- .fire_archive_manifest(zip)
  staging <- file.path(dest, sprintf(".staging-%d-%s", Sys.getpid(), basename(tempfile(""))))
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  .fire_archive_extract(zip, staging)
  short <- .fire_archive_shortfall(staging, manifest)
  if (length(short)) {
    stop(
      sprintf(
        "extraction of %s is incomplete: %d file(s) missing or short (e.g. %s)",
        basename(url),
        length(short),
        paste(utils::head(short, 3L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  .fire_archive_publish(staging, dest)
  try(
    writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), .fire_archive_stamp(dest, url)),
    silent = TRUE
  )

  found <- .fire_archive_files(dest, pattern)
  if (!length(found)) {
    stop(
      sprintf("no file matching '%s' after extracting %s", pattern, basename(url)),
      call. = FALSE
    )
  }
  found
}

#' Download + load NFDB fire points, harmonised + clipped to a study area
#'
#' Downloads the Canadian National Fire Database (NFDB) fire-point archive from the CWFIS open-data
#' server (~42 MB), extracts it, and loads it via [load_nfdb_points()] (same `YEAR` + `SIZE_HA`
#' harmonisation and study-area clipping, with the NFDB `CAUSE` column and other source attributes
#' preserved).
#'
#' The archive is downloaded + extracted once per `dest` and reused afterwards, so repeated runs --
#' and concurrent workers sharing a `dest` on shared storage, which coordinate via a lock -- do not
#' re-fetch it. A cached extraction is reused only once it has been verified complete against the
#' archive's own manifest, since a download or extraction interrupted partway leaves plausible-looking
#' short files behind that a reader will happily accept.
#'
#' @param study_area Study area defining the output CRS + crop extent: a file path (vector or raster),
#'   `sf`, `SpatVector`, or `SpatRaster` (e.g. a simulation `flammableMap`). It must carry a CRS. A
#'   **vector** study area selects records by its geometry; a **`SpatRaster`** selects by its extent,
#'   and its `NA` cells do not narrow that selection -- so an irregular study area passed as a raster
#'   also returns records from the bounding box around it.
#' @param fire_years Integer vector of fire years to keep; `NULL` (default) keeps all years.
#' @param min_size_ha Minimum fire size in hectares to keep (default `1`); records with a smaller
#'   reported `SIZE_HA` are dropped, while records with a missing (`NA`) size are always kept. Pass
#'   `0` to retain all fires (e.g. when a downstream model filters by size itself).
#' @param dest Directory to download + extract the archive into (created if needed). Defaults to a
#'   session tempdir; point it at a persistent (ideally shared) location to cache across runs.
#' @param url URL of the NFDB point shapefile archive (zip). `NULL` (default) uses the CWFIS
#'   current-version archive.
#'
#' @returns A `SpatVector` of NFDB fire points cropped to `study_area`, as [load_nfdb_points()].
#'
#' @seealso [load_nfdb_points()]
#' @family fire-record loaders
#' @inheritSection load_nbac_polys Spatial pre-filtering
#' @export
fetch_nfdb_points <- function(
  study_area,
  fire_years = NULL,
  min_size_ha = 1,
  dest = file.path(tempdir(), "NFDB_point"),
  url = NULL
) {
  if (is.null(url)) {
    url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point_shp.zip"
  }
  shp <- .download_fire_archive(url, dest, pattern = "NFDB_point.*\\.shp$")
  load_nfdb_points(shp, study_area = study_area, fire_years = fire_years, min_size_ha = min_size_ha)
}

#' Download + load NFDB fire polygons, harmonised + clipped to a study area
#'
#' Downloads the Canadian National Fire Database (NFDB) fire-polygon archive from the CWFIS open-data
#' server (~780 MB), extracts it, and loads it via [load_nfdb_polys()]. The record ships as several
#' multi-year partitions with differing columns; all of them are passed to the loader together, which
#' binds them. Caching, verification and concurrent-worker behaviour are as for [fetch_nfdb_points()].
#'
#' Prefer NBAC ([fetch_nbac_polys()]); use NFDB polygons only to backfill years NBAC does not cover.
#'
#' @inheritParams fetch_nfdb_points
#' @param url URL of the NFDB polygon shapefile archive (zip). `NULL` (default) uses the CWFIS
#'   current-version archive.
#'
#' @returns A `SpatVector` of NFDB polygons cropped to `study_area`, as [load_nfdb_polys()].
#'
#' @seealso [load_nfdb_polys()]
#' @family fire-record loaders
#' @inheritSection load_nbac_polys Spatial pre-filtering
#' @export
fetch_nfdb_polys <- function(
  study_area,
  fire_years = NULL,
  min_size_ha = 1,
  dest = file.path(tempdir(), "NFDB_poly"),
  url = NULL
) {
  if (is.null(url)) {
    url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  }
  shp <- .download_fire_archive(url, dest, pattern = "NFDB_poly_.*\\.shp$")
  load_nfdb_polys(shp, study_area = study_area, fire_years = fire_years, min_size_ha = min_size_ha)
}

#' Download + load NBAC fire perimeters, harmonised + clipped to a study area
#'
#' Downloads the National Burned Area Composite (NBAC) archive from the CWFIS open-data server
#' (~1.2 GB), extracts it, and loads it via [load_nbac_polys()]. Caching, verification and
#' concurrent-worker behaviour are as for [fetch_nfdb_points()].
#'
#' Unlike the NFDB archives, which live at a mutable `current_version/` path, NBAC releases carry a
#' date in the filename. The default `url` therefore pins a specific release rather than tracking the
#' newest one, so a pipeline re-run fetches the same data; pass `url` to move to a newer release
#' deliberately (the available releases are listed at
#' <https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/>).
#'
#' @inheritParams fetch_nfdb_points
#' @param url URL of the NBAC shapefile archive (zip). `NULL` (default) uses the pinned release
#'   described above.
#'
#' @returns A `SpatVector` of NBAC perimeters cropped to `study_area`, as [load_nbac_polys()].
#'
#' @seealso [load_nbac_polys()]
#' @family fire-record loaders
#' @inheritSection load_nbac_polys Spatial pre-filtering
#' @export
fetch_nbac_polys <- function(
  study_area,
  fire_years = NULL,
  min_size_ha = 1,
  dest = file.path(tempdir(), "NBAC"),
  url = NULL
) {
  if (is.null(url)) {
    url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/NBAC_1972to2025_20260513_shp.zip"
  }
  shp <- .download_fire_archive(url, dest, pattern = "NBAC_.*\\.shp$")
  load_nbac_polys(shp, study_area = study_area, fire_years = fire_years, min_size_ha = min_size_ha)
}
