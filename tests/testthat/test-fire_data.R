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
