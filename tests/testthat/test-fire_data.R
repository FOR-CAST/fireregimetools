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
