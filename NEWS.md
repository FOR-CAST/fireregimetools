# fireregimetools 0.1.0.9000

- `load_nbac_polys()` and `load_nfdb_polys()` now repair invalid fire-perimeter geometries (via `spatialutils::repair_geoms()`, which passes only the invalid subset to `terra::makeValid()`) instead of dropping them, recovering perimeters (e.g. ~767 of NFDB's ~41k national polygons) that were previously discarded.

# fireregimetools 0.0.1

- Initial release. Arrow-native, memory-bounded fire-regime summaries, consolidating the reusable fire-side post-processing logic from the `burnSummaries` module and project pipelines (companion to nrvtools).
- `write_burn_parquet()` and `open_burn_dataset()` write and lazily read per-replicate fire tables as partitioned parquet (`replicate=<rep>/part-0.parquet`), published atomically so concurrent writers on an NFS mount never collide.
- `summarize_fire_sizes()` and `summarize_annual_area()` reduce across replicates by pushing the count + five-number summary (`min`/`q25`/`median`/`q75`/`max`, plus `mean`/`sd`/`se`/`ci`) down to Arrow compute, so the per-record rows are never all held in memory at once.
- `fire_cycle()` computes the landscape fire cycle (years) from a per-replicate annual-area-burned table.
- `fire_size_histogram()` draws the fire-size distribution (count histogram with a median-log-size-per-bin overlay on a secondary axis).
