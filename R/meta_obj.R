#' MetaObject reference class
#'
#' A lightweight `setRefClass` container for storing a file path and a metadata
#' table extracted from that file/directory.
#'
#' When initialized with a non-empty `file_path`, the object stores the path in
#' `filepath` and populates `metadata` using \link{cfcore_extract_info}.
#'
#' @field filepath Character scalar. Path to the source directory (or file path).
#' @field metadata `data.frame`. Parsed metadata derived from `filepath`.
#'
#' @details
#' This is an R reference class (mutable). Instances are created via:
#' `MetaObject$new("path/to/dir")`.
#'
#' @seealso \link{cfcore_extract_info}
#' @keywords internal
mo <- setRefClass(
  "MetaObject",
  fields = list(
    filepath = "character",
    metadata = "data.frame"
  ),
  methods = list(
    initialize = function(file_path = character(0)) {
      if (base::length(file_path) > 0) {
        if (!is.character(file_path) || base::length(file_path) != 1L) {
          stop("`file_path` must be a single character path.", call. = FALSE)
        }

        .self$filepath <- file_path
        .self$metadata <- cfcore_extract_info(file_path)
      }
    }
  )
)


