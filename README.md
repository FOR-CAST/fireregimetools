# fireregimetools

<!-- badges: start -->
[![R-CMD-check](https://github.com/FOR-CAST/fireregimetools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FOR-CAST/fireregimetools/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Arrow-native, memory-bounded summaries of simulated fire regimes.

`fireregimetools` is the fire-side companion to
[`nrvtools`](https://github.com/FOR-CAST/nrvtools): it consolidates the reusable
fire-regime post-processing logic (fire-size distributions, annual area burned,
fire cycle) into a shared, SpaDES-agnostic package so it is implemented once
rather than reimplemented in each project.

Each replicate's fire-size or annual-area-burned table is written to its own
partition of an on-disk Arrow dataset (`write_burn_parquet()`), and the
across-replicate reduction is computed by pushing the aggregation down to Arrow
compute (`summarize_fire_sizes()`, `summarize_annual_area()`, `fire_cycle()`), so
the per-record rows are never all held in memory at once. `fire_size_histogram()`
draws the fire-size distribution.

## Installation

```r
renv::install("FOR-CAST/fireregimetools")
```

## License

Apache License (>= 2). See `LICENSE.md`.
