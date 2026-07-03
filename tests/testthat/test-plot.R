test_that("fire_size_histogram returns a ggplot from a vector or data frame", {
  withr::local_seed(1)
  sizes <- stats::rlnorm(200, meanlog = 3, sdlog = 1.5)

  p_vec <- fire_size_histogram(sizes)
  expect_s3_class(p_vec, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_vec))

  p_df <- fire_size_histogram(data.frame(size = sizes), size_col = "size")
  expect_s3_class(p_df, "ggplot")
})

test_that("fire_size_histogram returns NULL for empty or non-positive input", {
  expect_null(fire_size_histogram(numeric(0)))
  expect_null(fire_size_histogram(c(0, -1, NA)))
  expect_null(fire_size_histogram(data.frame(size = numeric(0))))
})
