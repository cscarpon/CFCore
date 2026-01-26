#' Extract file metadata from a directory
#'
#' Scans a directory and returns a data frame describing files with supported
#' extensions (LiDAR, rasters, vectors, tabular, and common serialized formats).
#'
#' @param file_path Character scalar. Directory path to scan.
#' @param exts_regex Regular expression used to match file extensions.
#'   Defaults to common CloudFlux inputs.
#'
#' @return A `data.frame` with columns:
#' - `id` integer row id
#' - `file_path` full normalized path
#' - `file_name` base name
#' - `size_mb` numeric size in MB
#' - `ext` lower-case file extension
#' - `mtime` file modification time (`POSIXct`)
#'
#' @keywords internal
cfcore_extract_info <- function(
    file_path,
    exts_regex = "\\.(laz|las|xyz|csv|shp|tif|tiff|json|rds)$"
) {
  if (!is.character(file_path) || base::length(file_path) != 1L) {
    stop("`file_path` must be a single character path.", call. = FALSE)
  }

  path_normal <- base::normalizePath(file_path, mustWork = FALSE)

  if (!base::file.exists(path_normal)) {
    stop("Path does not exist: ", path_normal, call. = FALSE)
  }
  if (!base::dir.exists(path_normal)) {
    stop("`file_path` must be a directory: ", path_normal, call. = FALSE)
  }

  data_list <- base::list.files(path_normal, pattern = exts_regex, full.names = TRUE)

  if (base::length(data_list) == 0L) {
    return(base::data.frame(
      id = integer(0),
      file_path = character(0),
      file_name = character(0),
      size_mb = numeric(0),
      ext = character(0),
      mtime = as.POSIXct(character(0)),
      stringsAsFactors = FALSE
    ))
  }

  rows <- base::lapply(base::seq_along(data_list), function(i) {
    current_file_path <- base::normalizePath(data_list[[i]], mustWork = FALSE)

    info <- base::file.info(current_file_path)
    ext <- base::tolower(tools::file_ext(current_file_path))

    base::data.frame(
      id = i,
      file_path = current_file_path,
      file_name = base::basename(current_file_path),
      size_mb = base::round(info$size / (1024^2), 2),
      ext = ext,
      mtime = info$mtime,
      stringsAsFactors = FALSE
    )
  })

  meta_df <- base::do.call(base::rbind, rows)

  meta_df <- meta_df |>
    dplyr::arrange(.data$ext) |>
    dplyr::mutate(id = dplyr::row_number()) |>
    dplyr::select(.data$id, dplyr::everything())

  meta_df
}

# Backward-compatible alias (optional).
# If you want to force callers to migrate, delete this alias.
#' @keywords internal
extract_info <- cfcore_extract_info

#' Capture console output from an expression
#'
#' Evaluates an expression and captures anything printed to stdout.
#'
#' @param expr An R expression to evaluate.
#' @return Character vector of captured lines.
#' @keywords internal
capture_output <- function(expr) {
  temp <- base::tempfile()
  base::sink(temp)
  on.exit(base::sink(), add = TRUE)
  on.exit(base::unlink(temp), add = TRUE)

  base::eval.parent(base::substitute(expr))

  base::readLines(temp, warn = FALSE)
}

#' Create a directory (recursively) if it does not exist
#'
#' @param in_dir Character scalar. Directory path to create.
#' @param verbose Logical. If `TRUE`, prints a message.
#' @return Invisibly returns `in_dir`.
#' @keywords internal
create_directories <- function(in_dir, verbose = FALSE) {
  if (!is.character(in_dir) || base::length(in_dir) != 1L) {
    stop("`in_dir` must be a single character path.", call. = FALSE)
  }

  if (!base::dir.exists(in_dir)) {
    base::dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(verbose)) base::message("Created directory: ", in_dir)
  } else {
    if (isTRUE(verbose)) base::message("Directory already exists: ", in_dir)
  }

  invisible(in_dir)
}

#' Time an expression
#'
#' Evaluates an expression and returns its result along with elapsed time.
#'
#' @param expr An R expression to evaluate.
#' @param verbose Logical. If `TRUE`, prints a message with elapsed seconds.
#' @return A list with `value` (the expression result) and `elapsed` (difftime).
#' @keywords internal
count_time <- function(expr, verbose = FALSE) {
  start_time <- base::Sys.time()
  value <- base::eval.parent(base::substitute(expr))
  end_time <- base::Sys.time()

  elapsed <- end_time - start_time

  if (isTRUE(verbose)) {
    base::message("Task took ", base::round(base::as.numeric(elapsed, units = "secs"), 2), " seconds.")
  }

  list(value = value, elapsed = elapsed)
}

#' Delete files in a directory
#'
#' Deletes files in a directory. Use with care.
#'
#' @param dir Character scalar. Directory to delete files from.
#' @param pattern Optional regex pattern to filter files. Defaults to `NULL` (all files).
#' @param recursive Logical. If `TRUE`, includes files in subdirectories.
#' @param dry_run Logical. If `TRUE`, returns the files that would be deleted without deleting.
#'
#' @return Character vector of files deleted (or that would be deleted if `dry_run=TRUE`).
#' @keywords internal
delete_all_files <- function(dir, pattern = NULL, recursive = FALSE, dry_run = FALSE) {
  if (!is.character(dir) || base::length(dir) != 1L) {
    stop("`dir` must be a single character path.", call. = FALSE)
  }
  if (!base::dir.exists(dir)) {
    return(character(0))
  }

  files <- base::list.files(dir, full.names = TRUE, recursive = recursive, pattern = pattern)

  if (base::length(files) == 0L) {
    return(character(0))
  }

  if (isTRUE(dry_run)) {
    return(files)
  }

  ok <- base::file.remove(files)
  files[ok]
}
