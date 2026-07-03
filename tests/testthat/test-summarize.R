test_that("summarize_fire_sizes matches an in-memory oracle and derives se", {
  root <- withr::local_tempdir()
  withr::local_seed(1)
  vals <- lapply(1:5, function(r) {
    df <- data.frame(scenario = "hrv", size = stats::runif(4, 1, 100))
    write_burn_parquet(df, root, r)
    df
  })
  env <- summarize_fire_sizes(root)
  expect_setequal(
    names(env),
    c("scenario", "n_fires", "mean", "sd", "min", "q25", "median", "q75", "max", "se", "ci")
  )
  ov <- do.call(rbind, vals)$size
  expect_equal(env$n_fires, length(ov))
  expect_equal(env$mean, mean(ov))
  expect_equal(env$min, min(ov))
  expect_equal(env$max, max(ov))
  expect_equal(env$se, env$sd / sqrt(env$n_fires))
})

test_that("summarize_annual_area groups per time and region and counts reps", {
  root <- withr::local_tempdir()
  withr::local_seed(2)
  for (r in 1:4) {
    write_burn_parquet(
      data.frame(
        scenario = "hrv",
        time = rep(c(0L, 100L), each = 2L),
        region = rep(c("all", "core"), times = 2L),
        area_ha = stats::runif(4, 0, 500)
      ),
      root,
      r
    )
  }
  env <- summarize_annual_area(root)
  expect_true(all(c("time", "region", "n_reps", "mean", "median") %in% names(env)))
  expect_equal(nrow(env), 4L) ## 2 times x 2 regions
  expect_equal(unique(env$n_reps), 4L)
})

test_that("fire_cycle = landscape_ha / mean annual area", {
  root <- withr::local_tempdir()
  for (r in 1:3) {
    write_burn_parquet(data.frame(scenario = "hrv", area_ha = c(100, 300)), root, r)
  }
  fc <- fire_cycle(root, landscape_ha = 6000)
  expect_equal(fc$mean_annual_area, 200)
  expect_equal(fc$fire_cycle_yr, 6000 / 200)
})

test_that("summarize_fire_sizes errors on a missing size column", {
  root <- withr::local_tempdir()
  write_burn_parquet(data.frame(scenario = "hrv", size = 1), root, 1L)
  expect_snapshot(summarize_fire_sizes(root, size_col = "nope"), error = TRUE)
})

test_that("empty input returns an empty data frame", {
  expect_equal(nrow(summarize_fire_sizes(character(0))), 0L)
  expect_equal(nrow(summarize_annual_area(character(0))), 0L)
  expect_equal(nrow(fire_cycle(character(0), landscape_ha = 1)), 0L)
})
