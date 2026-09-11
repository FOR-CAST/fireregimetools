## a small study-area polygon + a rasterToMatch over the same extent (EPSG:3005)
make_sa_vect <- function() {
  terra::vect("POLYGON ((0 0, 300 0, 300 300, 0 300, 0 0))", crs = "EPSG:3005")
}
make_sa_rast <- function() {
  terra::rast(make_sa_vect(), resolution = 30)
}

sq <- function(x0, y0, s = 100) {
  terra::vect(
    sprintf(
      "POLYGON ((%1$s %2$s, %3$s %2$s, %3$s %4$s, %1$s %4$s, %1$s %2$s))",
      x0,
      y0,
      x0 + s,
      y0 + s
    ),
    crs = "EPSG:3005"
  )
}

test_that("load_nbac_polys() tolerates alternate year/size columns + filters fire years/size", {
  nbac <- rbind(sq(0, 0), sq(120, 0), sq(0, 120))
  nbac$FIRE_YEAR <- c(2010L, 1999L, 2012L) # FIRE_YEAR (not YEAR)
  nbac$POLY_HA <- c(50, 200, 0.5) # POLY_HA (not ADJ_HA); 0.5 ha is below the >= 1 ha cutoff
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nbac, f, overwrite = TRUE)

  out <- load_nbac_polys(f, make_sa_vect(), fire_years = 2000:2020)
  ## 1999 out of range, 2012 below 1 ha -> only 2010 survives; harmonised to YEAR + SIZE_HA
  expect_s4_class(out, "SpatVector")
  expect_equal(out$YEAR, 2010L)
  expect_equal(out$SIZE_HA, 50)
})

test_that("study_area may be a SpatRaster (e.g. a flammableMap)", {
  nbac <- sq(0, 0)
  nbac$YEAR <- 2005L
  nbac$ADJ_HA <- 10
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nbac, f, overwrite = TRUE)

  out <- load_nbac_polys(f, make_sa_rast(), fire_years = 2000:2020)
  expect_s4_class(out, "SpatVector")
  expect_equal(out$YEAR, 2005L)
  expect_true(terra::same.crs(out, make_sa_rast()))
})

test_that("load_nfdb_polys() harmonises YEAR + SIZE_HA", {
  nfdb <- rbind(sq(0, 0), sq(120, 0))
  nfdb$YEAR <- c(1985L, 2015L)
  nfdb$SIZE_HA <- c(5, 8)
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  out <- load_nfdb_polys(f, make_sa_vect(), fire_years = 1900:2025)
  expect_setequal(out$YEAR, c(1985L, 2015L))
  expect_setequal(out$SIZE_HA, c(5, 8))
})

## a small NFDB-style point layer (x/y in EPSG:3005)
pts <- function(x, y, YEAR, SIZE_HA) {
  terra::vect(
    data.frame(x = x, y = y, YEAR = YEAR, SIZE_HA = SIZE_HA),
    geom = c("x", "y"),
    crs = "EPSG:3005"
  )
}

test_that("load_nfdb_points() harmonises YEAR + SIZE_HA and filters years/size", {
  p <- pts(
    x = c(50, 150, 250),
    y = c(50, 50, 50),
    YEAR = c(2010L, 1990L, 2011L),
    SIZE_HA = c(3, 12, 0.4) # 0.4 ha is below the >= 1 ha cutoff
  )
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(p, f, overwrite = TRUE)

  out <- load_nfdb_points(f, make_sa_vect(), fire_years = 2000:2020)
  ## 1990 out of range, 2011 below 1 ha -> only 2010 survives
  expect_s4_class(out, "SpatVector")
  expect_equal(terra::geomtype(out), "points")
  expect_equal(out$YEAR, 2010L)
  expect_equal(out$SIZE_HA, 3)
})

test_that("load_nfdb_points() min_size_ha controls the small-fire cutoff", {
  p <- pts(x = c(50, 150), y = c(50, 50), YEAR = c(2010L, 2011L), SIZE_HA = c(0.2, 3))
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(p, f, overwrite = TRUE)

  ## default 1 ha floor drops the 0.2 ha fire; min_size_ha = 0 keeps both
  expect_equal(load_nfdb_points(f, make_sa_vect(), fire_years = 2000:2020)$SIZE_HA, 3)
  expect_setequal(
    load_nfdb_points(f, make_sa_vect(), fire_years = 2000:2020, min_size_ha = 0)$SIZE_HA,
    c(0.2, 3)
  )
})

test_that("fire_years = NULL keeps all years (no year filter)", {
  p <- pts(
    x = c(50, 150, 250),
    y = c(50, 50, 50),
    YEAR = c(1985L, 2005L, 2020L),
    SIZE_HA = c(3, 4, 5)
  )
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(p, f, overwrite = TRUE)

  ## default fire_years = NULL -> no year filter; all three survive
  out <- load_nfdb_points(f, make_sa_vect())
  expect_setequal(out$YEAR, c(1985L, 2005L, 2020L))
})

test_that("load_nfdb_points() projects + crops to the study area", {
  p <- pts(
    x = c(50, 500), # 500 is outside the 0-300 study-area extent
    y = c(50, 50),
    YEAR = c(2005L, 2006L),
    SIZE_HA = c(4, 4)
  )
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(p, f, overwrite = TRUE)

  out <- load_nfdb_points(f, make_sa_rast(), fire_years = 2000:2020)
  expect_equal(out$YEAR, 2005L)
  expect_true(terra::same.crs(out, make_sa_rast()))
})

test_that("load_nbac_polys() errors when year/size columns are absent", {
  bad <- sq(0, 0)
  bad$SOMETHING <- 1L
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(bad, f, overwrite = TRUE)

  expect_error(
    load_nbac_polys(f, make_sa_vect(), fire_years = 2000:2020),
    "missing expected year/size columns"
  )
})

test_that("fetch_nfdb_points() reuses an already-extracted archive (no download)", {
  p <- pts(x = c(50, 250), y = c(50, 50), YEAR = c(2010L, 1990L), SIZE_HA = c(3, 4))
  dest <- withr::local_tempdir()
  terra::writeVector(p, file.path(dest, "NFDB_point_20250101.shp"), overwrite = TRUE)

  ## the .shp is already in `dest`, so fetch must load it WITHOUT touching `url` (which is unreachable)
  out <- fetch_nfdb_points(
    make_sa_vect(),
    fire_years = 2000:2020,
    dest = dest,
    url = "https://example.invalid/NFDB_point_shp.zip"
  )
  expect_s4_class(out, "SpatVector")
  expect_equal(terra::geomtype(out), "points")
  expect_equal(out$YEAR, 2010L) # 1990 filtered out
  expect_equal(out$SIZE_HA, 3)
})

test_that("fetch_nfdb_points() extracts a cached zip archive", {
  skip_if(!nzchar(Sys.which("zip")), "no `zip` binary")
  p <- pts(x = c(50, 250), y = c(50, 50), YEAR = c(2010L, 1990L), SIZE_HA = c(3, 4))
  shpdir <- withr::local_tempdir()
  terra::writeVector(p, file.path(shpdir, "NFDB_point_20250101.shp"), overwrite = TRUE)

  dest <- withr::local_tempdir()
  url <- "https://example.invalid/NFDB_point_shp.zip"
  withr::with_dir(shpdir, utils::zip(file.path(dest, basename(url)), list.files(), flags = "-q"))

  ## no .shp in `dest`, but the zip is -> fetch unzips it (no download) then loads
  out <- fetch_nfdb_points(make_sa_vect(), fire_years = 2000:2020, dest = dest, url = url)
  expect_s4_class(out, "SpatVector")
  expect_equal(out$YEAR, 2010L)
})

## a file:// URL for a local path, on every platform: POSIX paths are already absolute
## ("/tmp/x" -> "file:///tmp/x"), while a Windows drive letter needs the extra leading slash
## ("C:/Temp/x" -> "file:///C:/Temp/x")
file_url <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(p, "/")) {
    p <- paste0("/", p)
  }
  paste0("file://", p)
}

## zip the contents of `dir` (flat, as the CWFIS archives are) into the absolute path `zipfile`
zip_dir <- function(dir, zipfile) {
  withr::with_dir(dir, utils::zip(zipfile, list.files(), flags = "-q"))
  zipfile
}

## two NFDB-style polygon partitions, as the poly record actually ships
write_poly_partitions <- function(dir) {
  a <- sq(0, 0)
  a$YEAR <- 1985L
  a$SIZE_HA <- 5
  b <- sq(120, 0)
  b$YEAR <- 2022L
  b$SIZE_HA <- 8
  terra::writeVector(a, file.path(dir, "NFDB_poly_1972to2020_20240101.shp"), overwrite = TRUE)
  terra::writeVector(b, file.path(dir, "NFDB_poly_2021to2024_20250101.shp"), overwrite = TRUE)
  invisible(dir)
}

test_that("fetch_nfdb_polys() binds every partition in the archive, not just the first", {
  skip_if(!nzchar(Sys.which("zip")), "no `zip` binary")
  src <- withr::local_tempdir()
  write_poly_partitions(src)
  dest <- withr::local_tempdir()
  url <- "https://example.invalid/NFDB_poly.zip"
  zip_dir(src, file.path(dest, basename(url)))

  out <- fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  expect_setequal(out$YEAR, c(1985L, 2022L))
  expect_setequal(out$SIZE_HA, c(5, 8))
})

test_that("a verified extraction is stamped and reused", {
  skip_if(!nzchar(Sys.which("zip")), "no `zip` binary")
  src <- withr::local_tempdir()
  write_poly_partitions(src)
  dest <- withr::local_tempdir()
  url <- "https://example.invalid/NFDB_poly.zip"
  zip_dir(src, file.path(dest, basename(url)))

  fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  stamp <- file.path(dest, ".NFDB_poly.zip.complete")
  expect_true(file.exists(stamp))

  ## the stamp short-circuits verification: an unreadable archive is never consulted again
  writeLines("not a zip", file.path(dest, basename(url)))
  out <- fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  expect_setequal(out$YEAR, c(1985L, 2022L))
})

test_that("a truncated extraction is re-extracted, not accepted", {
  skip_if(!nzchar(Sys.which("zip")), "no `zip` binary")
  src <- withr::local_tempdir()
  write_poly_partitions(src)
  dest <- withr::local_tempdir()
  url <- "https://example.invalid/NFDB_poly.zip"
  zip_dir(src, file.path(dest, basename(url)))

  fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  truncated <- file.path(dest, "NFDB_poly_2021to2024_20250101.dbf")
  full_size <- file.size(truncated)
  writeBin(raw(8), truncated)
  unlink(file.path(dest, ".NFDB_poly.zip.complete")) # as if the earlier run had been interrupted

  out <- fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  expect_equal(file.size(truncated), full_size)
  expect_setequal(out$YEAR, c(1985L, 2022L))
})

test_that("a truncated archive is re-downloaded rather than extracted", {
  skip_if(!nzchar(Sys.which("zip")), "no `zip` binary")
  src <- withr::local_tempdir()
  write_poly_partitions(src)
  upstream <- withr::local_tempdir()
  good_zip <- zip_dir(src, file.path(upstream, "NFDB_poly.zip"))
  url <- file_url(good_zip)

  dest <- withr::local_tempdir()
  ## a partial download from an earlier run: right name, unreadable central directory
  writeBin(readBin(good_zip, "raw", 512L), file.path(dest, "NFDB_poly.zip"))

  out <- fetch_nfdb_polys(make_sa_vect(), dest = dest, url = url)
  expect_setequal(out$YEAR, c(1985L, 2022L))
  expect_equal(file.size(file.path(dest, "NFDB_poly.zip")), file.size(good_zip))
})

test_that("fetch_nbac_polys() reuses an already-extracted archive (no download)", {
  nbac <- sq(0, 0)
  nbac$YEAR <- 2005L
  nbac$ADJ_HA <- 10
  dest <- withr::local_tempdir()
  terra::writeVector(nbac, file.path(dest, "NBAC_1972to2025_20260513.shp"), overwrite = TRUE)

  out <- fetch_nbac_polys(
    make_sa_vect(),
    dest = dest,
    url = "https://example.invalid/NBAC_1972to2025_20260513_shp.zip"
  )
  expect_equal(out$YEAR, 2005L)
  expect_equal(out$SIZE_HA, 10)
})

test_that("a study area without a CRS errors instead of silently returning nothing", {
  nfdb <- sq(0, 0)
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 5
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  ## note terra::rast() silently assigns WGS84 when the extent fits inside lon/lat bounds, so this
  ## needs an extent that does not
  nocrs <- terra::rast(terra::ext(0, 300, 0, 300), resolution = 30)
  expect_error(load_nfdb_polys(f, nocrs), "no CRS")
})

test_that("records without a CRS name the offending file", {
  nfdb <- sq(0, 0)
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 5
  ## a shapefile with its .prj removed, NOT a CRS-less GPKG: GPKG always records an SRS entry, and
  ## what an undefined one reads back as varies with the GDAL/PROJ build
  dir <- withr::local_tempdir()
  f <- file.path(dir, "NFDB_poly_nocrs.shp")
  terra::writeVector(nfdb, f, overwrite = TRUE)
  unlink(file.path(dir, "NFDB_poly_nocrs.prj"))
  expect_equal(terra::crs(terra::vect(f)), "") # the fixture is only useful if this holds

  expect_error(load_nfdb_polys(f, make_sa_vect()), "fire records have no CRS")
})

test_that("an empty result keeps the harmonised schema", {
  nfdb <- sq(0, 0)
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 5
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  ## terra::crop() hands back an attribute-less SpatVector when nothing survives; callers reaching
  ## for `$YEAR` would get NULL rather than an empty vector
  disjoint <- terra::vect(
    "POLYGON ((1e6 1e6, 2e6 1e6, 2e6 2e6, 1e6 2e6, 1e6 1e6))",
    crs = "EPSG:3005"
  )
  out <- load_nfdb_polys(f, disjoint)
  expect_equal(nrow(out), 0L)
  expect_true(all(c("YEAR", "SIZE_HA") %in% names(out)))
  expect_length(out$YEAR, 0L)
})

test_that("a vector study area selects by geometry; a SpatRaster selects by extent", {
  ## an L-shaped study area: the notch is inside the bounding box but outside the polygon
  sa <- terra::vect(
    "POLYGON ((0 0, 300 0, 300 150, 150 150, 150 300, 0 300, 0 0))",
    crs = "EPSG:3005"
  )
  nfdb <- rbind(sq(10, 10, s = 50), sq(200, 200, s = 50)) # inside; in the notch
  nfdb$YEAR <- c(2010L, 2011L)
  nfdb$SIZE_HA <- c(5, 5)
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  expect_equal(load_nfdb_polys(f, sa)$YEAR, 2010L)
  expect_setequal(load_nfdb_polys(f, terra::rast(sa, resolution = 30))$YEAR, c(2010L, 2011L))
})

test_that("a clipped edge perimeter keeps its full reported SIZE_HA", {
  nfdb <- sq(250, 100) # 100 x 100, straddling the eastern edge of the 0-300 study area
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 1000 # reported size of the whole fire
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  out <- load_nfdb_polys(f, make_sa_vect())
  expect_equal(out$SIZE_HA, 1000) # attribute untouched ...
  expect_equal(terra::expanse(out), 100 * 50) # ... while half the geometry was clipped away
})

## --- spatial pre-filtering (pushdown) ----------------------------------------------------------
##
## The pushdown is an optimisation, so what these tests pin down is that it changes nothing: every
## loader is run with `fireregimetools.prefilter` on and off and the two results compared. The
## `FALSE` path is the pre-pushdown code (each source read whole), so it is the "before" to the
## pushdown's "after".

## everything about a loaded record set that a caller can observe
fire_digest <- function(v) {
  list(
    n = nrow(v),
    names = sort(names(v)),
    crs = terra::crs(v),
    geom = if (nrow(v) > 0L) terra::geom(v) else NULL,
    years = v$YEAR,
    sizes = v$SIZE_HA
  )
}

## the same loader call with the pushdown on and off
both_ways <- function(fun, ...) {
  off <- withr::with_options(list(fireregimetools.prefilter = FALSE), fun(...))
  on <- withr::with_options(list(fireregimetools.prefilter = TRUE), fun(...))
  list(off = fire_digest(off), on = fire_digest(on))
}

## fire polygons on a lon/lat grid across Canada, written once per test that needs them: a lon/lat
## source read against a projected study area is the case that stresses the pushdown hardest, since
## the study area's boundary is a curve in the source CRS
lonlat_fire_grid <- function(path, nx = 24L, ny = 12L) {
  xs <- seq(-135, -60, length.out = nx)
  ys <- seq(48, 68, length.out = ny)
  cells <- expand.grid(x = xs, y = ys)
  polys <- lapply(seq_len(nrow(cells)), function(i) {
    terra::vect(
      sprintf(
        "POLYGON ((%1$f %2$f, %3$f %2$f, %3$f %4$f, %1$f %4$f, %1$f %2$f))",
        cells$x[i],
        cells$y[i],
        cells$x[i] + 1.2,
        cells$y[i] + 0.9
      ),
      crs = "EPSG:4326"
    )
  })
  v <- do.call(rbind, polys)
  v$YEAR <- as.integer(1972 + seq_len(nrow(v)) %% 50L)
  v$POLY_HA <- 10 * seq_len(nrow(v)) # an area column both the NBAC and NFDB loaders accept
  terra::writeVector(v, path, overwrite = TRUE)
  v
}

test_that("the pushdown selects a subset of the source, not all of it", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  src <- lonlat_fire_grid(f)
  sa <- terra::project(
    terra::vect("POLYGON ((-116 53, -112 53, -112 56, -116 56, -116 53))", crs = "EPSG:4326"),
    "EPSG:3978"
  )

  ## the pushdown must actually be pushing something down: a study area covering a few grid cells
  ## has to read far fewer than the source's features, or the optimisation has silently regressed
  read <- .read_fire_source(f, sa)
  expect_lt(nrow(read), nrow(src) / 4)
  expect_gt(nrow(read), 0L)
})

test_that("the pushdown returns the same records as reading the sources whole", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  lonlat_fire_grid(f)
  sa <- terra::project(
    terra::vect("POLYGON ((-116 53, -112 53, -112 56, -116 56, -116 53))", crs = "EPSG:4326"),
    "EPSG:3978"
  )

  r <- both_ways(load_nbac_polys, f, sa, fire_years = NULL, min_size_ha = 0)
  expect_gt(r$on$n, 0L)
  expect_equal(r$on, r$off)
})

test_that("the pushdown is invariant across randomised study areas", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  lonlat_fire_grid(f)
  withr::local_seed(20260911)

  ## random study areas of random sizes, in a CRS the source is not in: each must return exactly
  ## what a whole-source read would, including the ones that select nothing
  selected <- integer(0)
  for (i in seq_len(25L)) {
    x0 <- stats::runif(1, -134, -62)
    y0 <- stats::runif(1, 48, 67)
    w <- stats::runif(1, 0.2, 8)
    h <- stats::runif(1, 0.2, 5)
    sa <- terra::project(
      terra::vect(
        sprintf(
          "POLYGON ((%1$f %2$f, %3$f %2$f, %3$f %4$f, %1$f %4$f, %1$f %2$f))",
          x0,
          y0,
          x0 + w,
          y0 + h
        ),
        crs = "EPSG:4326"
      ),
      "EPSG:3978"
    )
    r <- both_ways(load_nfdb_polys, f, sa, min_size_ha = 0)
    expect_equal(r$on, r$off, info = sprintf("study area %d: %f %f %f %f", i, x0, y0, w, h))
    selected <- c(selected, r$on$n)
  }

  ## an agreement over nothing but empty results would prove nothing
  expect_gt(sum(selected > 0L), 15L)
  expect_gt(max(selected), 1L)
})

test_that("the pushdown is invariant for a SpatRaster study area (extent semantics)", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  lonlat_fire_grid(f)
  ## an L-shaped footprint: as a raster it selects by extent, so records in the notch are kept --
  ## the pushdown must not narrow that to the polygon
  sa <- terra::project(
    terra::vect(
      "POLYGON ((-124 50, -100 50, -100 56, -112 56, -112 62, -124 62, -124 50))",
      crs = "EPSG:4326"
    ),
    "EPSG:3978"
  )
  r_vect <- both_ways(load_nfdb_polys, f, sa, min_size_ha = 0)
  r_rast <- both_ways(load_nfdb_polys, f, terra::rast(sa, resolution = 5000), min_size_ha = 0)

  expect_equal(r_vect$on, r_vect$off)
  expect_equal(r_rast$on, r_rast$off)
  expect_gt(r_rast$on$n, r_vect$on$n) # the notch really is populated
})

test_that("the pushdown keeps a perimeter straddling the study-area edge", {
  ## the pushdown's one failure mode would be dropping a record that the crop keeps, and an edge
  ## record read in a different CRS is where that would happen
  nfdb <- terra::vect(
    "POLYGON ((-112.5 53.5, -110 53.5, -110 55, -112.5 55, -112.5 53.5))",
    crs = "EPSG:4326"
  )
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 50000
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  sa <- terra::project(
    terra::vect("POLYGON ((-116 53, -112 53, -112 56, -116 56, -116 53))", crs = "EPSG:4326"),
    "EPSG:3978"
  )
  r <- both_ways(load_nfdb_polys, f, sa)
  expect_equal(r$on$n, 1L)
  expect_equal(r$on, r$off)
})

test_that("the pushdown is invariant for point records", {
  p <- terra::vect(
    data.frame(
      x = seq(-130, -70, length.out = 60),
      y = seq(50, 65, length.out = 60),
      YEAR = as.integer(1972 + seq_len(60L) %% 50L),
      SIZE_HA = seq_len(60L)
    ),
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(p, f, overwrite = TRUE)

  sa <- terra::project(
    terra::vect("POLYGON ((-110 56, -95 56, -95 61, -110 61, -110 56))", crs = "EPSG:4326"),
    "EPSG:3978"
  )
  r <- both_ways(load_nfdb_points, f, sa, min_size_ha = 0)
  expect_gt(r$on$n, 0L)
  expect_equal(r$on, r$off)
})

test_that("a partition with no records near the study area binds with one that has them", {
  ## the NFDB polygon record ships multi-year partitions, and a study area will often sit inside
  ## only some of them; the empty reads must not cost the schema (or the records)
  dir <- withr::local_tempdir()
  near <- sq(0, 0)
  near$YEAR <- 1985L
  near$SIZE_HA <- 5
  far <- terra::vect(
    "POLYGON ((1e6 1e6, 1000100 1e6, 1000100 1000100, 1e6 1000100, 1e6 1e6))",
    crs = "EPSG:3005"
  )
  far$YEAR <- 2022L
  far$SIZE_HA <- 8
  fa <- file.path(dir, "NFDB_poly_a.gpkg")
  fb <- file.path(dir, "NFDB_poly_b.gpkg")
  terra::writeVector(near, fa, overwrite = TRUE)
  terra::writeVector(far, fb, overwrite = TRUE)

  r <- both_ways(load_nfdb_polys, c(fa, fb), make_sa_vect())
  expect_equal(r$on$n, 1L)
  expect_equal(r$on$years, 1985L)
  expect_equal(r$on, r$off)
})

test_that("a source with no CRS is reported before its geometry is read", {
  nfdb <- sq(0, 0)
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 5
  dir <- withr::local_tempdir()
  f <- file.path(dir, "NFDB_poly_nocrs.shp")
  terra::writeVector(nfdb, f, overwrite = TRUE)
  unlink(file.path(dir, "NFDB_poly_nocrs.prj"))

  ## the proxy carries the (missing) CRS, so the pushdown path raises this too -- and names the file
  expect_error(load_nfdb_polys(f, make_sa_vect()), "NFDB_poly_nocrs\\.shp")
})

test_that(".prefilter_extent() bounds the study area in the source CRS", {
  sa <- terra::project(
    terra::vect("POLYGON ((-116 53, -112 53, -112 56, -116 56, -116 53))", crs = "EPSG:4326"),
    "EPSG:3978"
  )
  e <- .prefilter_extent(sa, "EPSG:4326")

  ## it must cover the study area's true lon/lat footprint ...
  true_ext <- terra::ext(terra::project(sa, "EPSG:4326"))
  expect_equal(as.vector(terra::intersect(e, true_ext)), as.vector(true_ext))
  ## ... and still be a small part of a national source's extent
  expect_lt((e$xmax - e$xmin) * (e$ymax - e$ymin), 0.05 * (141 - 52) * (83 - 41))
})

test_that(".prefilter_extent() declines rather than guess", {
  sa <- terra::vect("POLYGON ((0 0, 300 0, 300 300, 0 300, 0 0))", crs = "EPSG:3005")

  expect_null(.prefilter_extent(sa, "")) # no source CRS to project into
  ## a degenerate (zero-area) study area
  expect_null(.prefilter_extent(terra::vect("POINT (10 10)", crs = "EPSG:3005"), "EPSG:4326"))
  ## a study area reaching outside the source CRS's domain: project() drops the points it cannot
  ## convert, so a box around the survivors would be too small -- decline instead
  world <- terra::vect(
    "POLYGON ((-179 -80, 179 -80, 179 80, -179 80, -179 -80))",
    crs = "EPSG:4326"
  )
  ortho <- "+proj=ortho +lat_0=60 +lon_0=-100"
  ## which holds only where the PROJ build refuses the far hemisphere, as they do by handing back
  ## empty geometries for those points -- the very thing the guard is there to catch
  pts <- terra::vect(
    as.matrix(expand.grid(x = seq(-179, 179, length.out = 8L), y = seq(-80, 80, length.out = 8L))),
    type = "points",
    crs = "EPSG:4326"
  )
  kept <- nrow(terra::crds(suppressWarnings(terra::project(pts, ortho))))
  skip_if(
    kept == nrow(terra::crds(pts)),
    "this PROJ build projects the whole world orthographically"
  )
  expect_null(.prefilter_extent(world, ortho))
})

test_that("prefilter = FALSE still reads (and errors) the same way", {
  withr::local_options(fireregimetools.prefilter = FALSE)
  nfdb <- sq(0, 0)
  nfdb$YEAR <- 2010L
  nfdb$SIZE_HA <- 5
  f <- withr::local_tempfile(fileext = ".gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  expect_equal(load_nfdb_polys(f, make_sa_vect())$YEAR, 2010L)

  nocrs <- terra::rast(terra::ext(0, 300, 0, 300), resolution = 30)
  expect_error(load_nfdb_polys(f, nocrs), "no CRS")
})

test_that("a source that cannot be pre-filtered is read whole, and says so", {
  nfdb <- rbind(sq(0, 0), sq(120, 0))
  nfdb$YEAR <- c(2010L, 2011L)
  nfdb$SIZE_HA <- c(5, 6)
  ## a fixed basename: it appears in the message below
  f <- file.path(withr::local_tempdir(), "NFDB_poly_x.gpkg")
  terra::writeVector(nfdb, f, overwrite = TRUE)

  ## a zero-area study area yields no usable filter extent, so the source is read whole -- the
  ## fallback is correct, but on the national records it costs minutes, so it is announced
  expect_snapshot(v <- .read_fire_source(f, terra::vect("POINT (10 10)", crs = "EPSG:3005")))
  expect_equal(nrow(v), 2L)
})
