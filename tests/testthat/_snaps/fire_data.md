# a source that cannot be pre-filtered is read whole, and says so

    Code
      v <- .read_fire_source(f, terra::vect("POINT (10 10)", crs = "EPSG:3005"))
    Message
      no spatial pre-filter for NFDB_poly_x.gpkg (the study area has no usable footprint in its CRS); reading the whole source

