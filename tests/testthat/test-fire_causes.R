## records with the cause and prescribed-burn columns of NBAC / NFDB, as the loaders keep them
fire_points <- function(n, ...) {
  v <- terra::vect(cbind(seq_len(n), seq_len(n)), crs = "EPSG:3978")
  cols <- list(...)
  for (nm in names(cols)) {
    v[[nm]] <- cols[[nm]]
  }
  v
}

test_that("NBAC causes fall in their groups, ignoring case; blank and missing are unknown", {
  d <- data.frame(FIRECAUS = c("Natural", "Human", "Undetermined", "", NA, " natural ", "Unknown"))
  expect_identical(
    fire_cause_group(d, "nbac"),
    c("natural", "human", "unknown", "unknown", "unknown", "natural", "unknown")
  )
})

test_that("NFDB causes fall in their groups, in the spellings of each release", {
  d <- data.frame(CAUSE = c("N", "L", "H", "H-PB", "U", "n/a", "Re", "RE", "", NA))
  expect_identical(
    fire_cause_group(d, "nfdb"),
    c(
      "natural",
      "natural",
      "human",
      "human",
      "unknown",
      "unknown",
      "other",
      "other",
      "unknown",
      "unknown"
    )
  )
})

test_that("a cause value not in the tables is an error naming it, not a guess", {
  expect_error(fire_cause_group(data.frame(FIRECAUS = c("Natural", "Arson")), "nbac"), "Arson")
  expect_error(fire_cause_group(data.frame(CAUSE = "X"), "nfdb"), "NFDB `CAUSE`.*X")
})

test_that("a record without its cause column is an error naming the column", {
  expect_error(fire_cause_group(data.frame(CAUSE = "N"), "nbac"), "FIRECAUS")
  expect_identical(
    fire_cause_group(data.frame(CAUSE2 = "N"), "nfdb", cause_col = "CAUSE2"),
    "natural"
  )
})

test_that("prescribed burns are found by PRESCRIBED, FIRE_TYPE, CAUSE and CAUSE2", {
  d <- data.frame(
    PRESCRIBED = c("true", "", NA, "No", "0", "PB", "1", "yes", "false", ""),
    FIRE_TYPE = c("", "", "", "", "", "", "", "", "", "PB"),
    CAUSE = c("N", "N", "H-PB", "H", "N", "H", "H", "H", "H", "H"),
    CAUSE2 = c("", "H-PB", "", "", "", "", "", "", "", "")
  )
  expect_identical(
    is_prescribed_burn(d),
    c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
})

test_that("columns the records lack are not consulted", {
  expect_identical(is_prescribed_burn(data.frame(YEAR = 2000:2001)), c(FALSE, FALSE))
  expect_identical(is_prescribed_burn(data.frame(PRESCRIBED = "true")), TRUE)
})

test_that("by default only prescribed burns are dropped, and no cause column is needed", {
  v <- fire_points(3, YEAR = 2000:2002, PRESCRIBED = c("", "true", ""))
  kept <- filter_fire_records(v)
  expect_s4_class(kept, "SpatVector")
  expect_identical(kept$YEAR, c(2000L, 2002L))
  expect_equal(nrow(filter_fire_records(v, drop_prescribed = FALSE)), 3)
})

test_that("causes keeps those groups and still drops prescribed burns", {
  v <- fire_points(
    5,
    YEAR = 2000:2004,
    FIRECAUS = c("Natural", "Human", "Undetermined", "Natural", "Human"),
    PRESCRIBED = c("", "", "", "true", "")
  )
  expect_identical(filter_fire_records(v, c("natural", "unknown"), "nbac")$YEAR, c(2000L, 2002L))
  expect_identical(filter_fire_records(v, "human", "nbac")$YEAR, c(2001L, 2004L))
  expect_identical(
    filter_fire_records(v, "natural", "nbac", drop_prescribed = FALSE)$YEAR,
    c(2000L, 2003L)
  )
})

test_that("NFDB records are selected by CAUSE", {
  d <- data.frame(YEAR = 1960:1964, CAUSE = c("N", "H", "H-PB", "U", "Re"))
  ## H-PB is a human cause, but a prescribed burn: dropped either way
  expect_identical(filter_fire_records(d, c("human", "other"), "nfdb")$YEAR, c(1961L, 1964L))
  expect_identical(filter_fire_records(d, source = "nfdb")$YEAR, c(1960L, 1961L, 1963L, 1964L))
})

test_that("an unknown cause group is an error listing the groups", {
  d <- data.frame(FIRECAUS = "Natural")
  expect_error(
    filter_fire_records(d, "lightning", "nbac"),
    "lightning.*natural, unknown, human, other"
  )
})

test_that("no record left keeps the columns", {
  v <- fire_points(2, YEAR = 2000:2001, FIRECAUS = c("Human", "Human"))
  kept <- filter_fire_records(v, "natural", "nbac")
  expect_equal(nrow(kept), 0)
  expect_true(all(c("YEAR", "FIRECAUS") %in% names(kept)))
})
