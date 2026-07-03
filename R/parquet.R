#' Write one replicate's fire table to a partitioned parquet
#'
#' Writes `df` (one replicate's fire-size or annual-area table) under
#' `root/replicate=<replicate>/part-0.parquet`, stamping `replicate` as a data
#' column. The write is atomic on POSIX/NFS: the parquet is written to a unique
#' temporary name in the destination directory and renamed into place, so a
#' concurrent reader (or a retried write) never observes a partial file. Many
#' replicate writers can therefore run concurrently against an NFS output
#' directory without contention (each writes its own `replicate=` partition).
#'
#' Run-level identifiers (`studyArea`, `scenario`, `simArea`, ...) should already
#' be columns of `df` so that the summary functions can group by them when a union
#' of roots is opened together.
#'
#' @param df A `data.frame` of one replicate's fire records. `NULL` or a zero-row
#'   input returns `NULL` (nothing is written).
#' @param root Directory under which the `replicate=<replicate>` partition is
#'   created (created recursively if needed).
#' @param replicate Replicate identifier (integer or string) used for the
#'   partition directory and stamped as the `replicate` column.
#'
#' @return The written parquet path (length-1 character), or `NULL` for empty
#'   input.
#'
#' @export
write_burn_parquet <- function(df, root, replicate) {
  if (is.null(df) || nrow(df) == 0L) {
    return(NULL)
  }
  df[["replicate"]] <- replicate
  dst_dir <- file.path(root, paste0("replicate=", replicate))
  fs::dir_create(dst_dir)
  dst <- file.path(dst_dir, "part-0.parquet")
  tmp <- tempfile(tmpdir = dst_dir, fileext = ".parquet.tmp")
  arrow::write_parquet(df, tmp)
  if (!file.rename(tmp, dst)) {
    unlink(tmp)
    stop("write_burn_parquet(): failed to publish ", dst, call. = FALSE)
  }
  dst
}

#' Open replicate fire parquet(s) as one lazy Arrow dataset
#'
#' Resolves `x` -- parquet file paths and/or dataset roots -- to a flat list of
#' `*.parquet` files and opens them as a single lazy [arrow::open_dataset()]
#' `Dataset`. Because `replicate` and any run-level identifiers are stored as data
#' columns (not inferred from the directory tree), files from several roots open
#' together without a `UnionDataset`, and opening from explicit file paths avoids
#' depending on NFS directory-listing freshness.
#'
#' @param x Character vector of parquet file paths and/or directories; directories
#'   are scanned recursively for `*.parquet`.
#'
#' @return An Arrow `Dataset` (lazy), or `NULL` if no parquet files are found.
#'
#' @export
open_burn_dataset <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) == 0L) {
    return(NULL)
  }
  is_dir <- dir.exists(x)
  files <- c(
    x[!is_dir & grepl("\\.parquet$", x)],
    if (any(is_dir)) {
      list.files(x[is_dir], pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
    }
  )
  files <- files[file.exists(files)]
  if (length(files) == 0L) {
    return(NULL)
  }
  arrow::open_dataset(files)
}
