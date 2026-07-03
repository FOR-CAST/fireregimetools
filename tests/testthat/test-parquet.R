test_that("write_burn_parquet writes an atomic hive partition and stamps replicate", {
  root <- withr::local_tempdir()
  p <- write_burn_parquet(data.frame(year = c(1, 2), size = c(10, 20)), root, 3L)
  expect_match(p, "replicate=3/part-0\\.parquet$")
  expect_length(list.files(root, pattern = "\\.tmp$", recursive = TRUE), 0L)
  d <- dplyr::collect(arrow::open_dataset(p))
  expect_equal(unique(d$replicate), 3L)
})

test_that("write_burn_parquet returns NULL for empty input", {
  root <- withr::local_tempdir()
  expect_null(write_burn_parquet(NULL, root, 1L))
  expect_null(write_burn_parquet(data.frame(), root, 1L))
})

test_that("open_burn_dataset opens files or roots equivalently, NULL when empty", {
  root <- withr::local_tempdir()
  for (r in 1:3) {
    write_burn_parquet(data.frame(year = 1, size = r), root, r)
  }
  files <- list.files(root, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  expect_equal(nrow(dplyr::collect(open_burn_dataset(root))), 3L)
  expect_equal(nrow(dplyr::collect(open_burn_dataset(files))), 3L)
  expect_null(open_burn_dataset(character(0)))
  expect_null(open_burn_dataset(withr::local_tempdir()))
})
