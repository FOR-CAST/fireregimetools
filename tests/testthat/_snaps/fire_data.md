# a source that cannot be pre-filtered is read whole, and says so

    Code
      v <- .read_fire_source(f, terra::vect("POINT (10 10)", crs = "EPSG:3005"))
    Message
      no spatial pre-filter for NFDB_poly_x.gpkg (the study area has no usable footprint in its CRS); reading the whole source

# include_pre1972 is rejected for records that have no pre-1972 split

    Code
      download_fire_records("nbac_poly", dest = tempdir(), include_pre1972 = TRUE)
    Condition
      Error:
      ! `include_pre1972` applies only to `record = "nfdb_poly"`; NBAC begins in 1972 and the NFDB point record is not split by year.

---

    Code
      download_fire_records("nfdb_point", dest = tempdir(), include_pre1972 = TRUE)
    Condition
      Error:
      ! `include_pre1972` applies only to `record = "nfdb_poly"`; NBAC begins in 1972 and the NFDB point record is not split by year.

