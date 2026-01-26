#' Spatial container for LiDAR point clouds, rasters, and vector layers
#'
#' A reference-class object used as a mutable state container for CloudFluxCore
#' workflows. The container can be initialized from a `.las` / `.laz` file or an
#' `.xyz` text file and stores derived products such as DTM and nDSM rasters,
#' spatial index geometry, and optional masks/buildings layers.
#'
#' This class is intended to hold intermediate objects during processing.
#' For API usage, prefer passing and returning plain objects (e.g., `LAS`,
#' `SpatRaster`, `sf`) from standalone functions, and use `spatial_container`
#' only as an orchestration/state convenience.
#'
#' @section Fields:
#' \describe{
#'   \item{xyz}{`data.frame`. Point data with columns `X`, `Y`, `Z` (if initialized from `.xyz`).}
#'   \item{LPC}{`lidR::LAS`. LiDAR point cloud.}
#'   \item{ndsm}{`terra::SpatRaster`. Derived nDSM, optionally clipped.}
#'   \item{ndsm_raw}{`terra::SpatRaster`. Unclipped or pre-processed nDSM intermediate.}
#'   \item{DTM}{`terra::SpatRaster`. Derived DTM, optionally clipped.}
#'   \item{DTM_raw}{`terra::SpatRaster`. Unclipped DTM intermediate.}
#'   \item{index}{`sf::sfc`. Spatial extent geometry (bounding box) for the point cloud.}
#'   \item{mask}{`sf`. Optional polygon mask used for clipping rasters.}
#'   \item{buildings}{`sf`. Optional building footprints.}
#'   \item{filepath}{`character`. Source file path provided at initialization.}
#'   \item{filename}{`character`. Basename of `filepath`.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{initialize(file_path = character(0))}{Initialize the container from a `.las`, `.laz`, or `.xyz` file.}
#'   \item{set_crs(crs)}{Assign/transform CRS on stored spatial objects.}
#'   \item{get_data()}{Return the `xyz` data.frame.}
#'   \item{get_lpc()}{Return the `LPC` LAS object.}
#'   \item{to_xyz(path)}{Write `xyz` to a space-delimited XYZ file.}
#'   \item{to_dtm(resolution = 1)}{Create a DTM from ground points and optionally clip by `mask`.}
#'   \item{to_ndsm(resolution = 1)}{Create an nDSM using the stored DTM.}
#'   \item{save_mask(path)}{Write `mask` to disk via `sf::st_write()`.}
#'   \item{save_las(path)}{Write `LPC` to disk via `lidR::writeLAS()`.}
#'   \item{save_dtm(path)}{Write `DTM` to disk via `terra::writeRaster()`.}
#'   \item{save_ndsm(path)}{Write `ndsm` to disk via `terra::writeRaster()`.}
#'   \item{save_sc(path)}{Serialize the full reference-class object via `save()`.}
#' }
#'
#' @name spatial_container
#' @docType class
#' @aliases spatial_container-class
#' @keywords internal
NULL

# Ensure required namespaces are loaded before defining RefClasses that reference their S4 classes
if (!base::requireNamespace("lidR", quietly = TRUE)) {
  stop("Package 'lidR' is required but not installed.", call. = FALSE)
}
if (!base::requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required but not installed.", call. = FALSE)
}

# Register sf S3 classes so RefClass fields can use them as class names
methods::setOldClass("sf")
methods::setOldClass("sfc")

#' @import methods
#' @importClassesFrom lidR LAS
#' @importClassesFrom terra SpatRaster
NULL

spatial_container <- setRefClass(
  "spatial_container",
  fields = list(
    xyz = "data.frame",
    LPC = "LAS",
    ndsm = "SpatRaster",
    ndsm_raw = "SpatRaster",
    DTM = "SpatRaster",
    DTM_raw = "SpatRaster",
    index = "sfc",
    mask = "sf",
    buildings = "sf",
    filepath = "character",
    filename = "character"
  ),
  methods = list(
    initialize = function(file_path = character(0)) {
      .self$filepath <- file_path
      .self$filename <- base::basename(file_path)

      ext <- tools::file_ext(file_path)

      dummy_spat <- terra::rast(
        xmin = 0, xmax = 1, ymin = 0, ymax = 1,
        resolution = 1, vals = NA
      )
      .self$DTM <- dummy_spat
      .self$ndsm <- dummy_spat

      if (length(file_path) == 0L || identical(file_path, character(0)) || is.na(file_path) || file_path == "") {
        return(invisible(NULL))
      }

      if (!base::file.exists(file_path)) {
        stop("File does not exist: ", file_path, call. = FALSE)
      }

      if (ext == "xyz") {
        xyz_table <- utils::read.table(file_path)
        base::names(xyz_table) <- c("X", "Y", "Z")
        .self$xyz <<- xyz_table

        las_xyz <- lidR::LAS(xyz_table)

        ground <- lidR::classify_ground(las_xyz, algorithm = lidR::pmf(ws = 5, th = 3))
        las_xyz@data$Classification <- ground@data$Classification
        las_xyz@data$ReturnNumber <- 1L
        las_xyz@data$NumberOfReturns <- 1L

        .self$LPC <<- las_xyz

        # bbox/index from the LAS we just created (was previously incorrectly using `las`)
        las_sf <- sf::st_as_sf(las_xyz@data, coords = c("X", "Y"))
        las_extent <- sf::st_as_sfc(sf::st_bbox(las_sf))

        # xyz has no CRS by default; keep NA unless you explicitly set later via set_crs()
        sf::st_crs(las_extent) <- sf::st_crs(las_sf)

        .self$index <<- las_extent
      } else if (ext == "las" || ext == "laz") {
        las <- lidR::readLAS(file_path)

        # lidR::readLAS can return empty/invalid LAS; fail early
        if (is.null(las) || !inherits(las, "LAS")) {
          stop("Failed to read LAS/LAZ file: ", file_path, call. = FALSE)
        }

        las_crs <- sf::st_crs(las)
        .self$LPC <<- las

        las_sf <- sf::st_as_sf(las@data, coords = c("X", "Y"), crs = las_crs)
        las_extent <- sf::st_as_sfc(sf::st_bbox(las_sf), crs = las_crs)

        .self$index <<- las_extent
      } else {
        stop("Unsupported file extension: .", ext, call. = FALSE)
      }

      invisible(NULL)
    },
    set_crs = function(crs) {
      crs <- base::as.integer(crs)

      current_crs <- sf::st_crs(.self$LPC)

      if (is.null(current_crs)) {
        base::cat("The LPC does not have an associated CRS.\n")
        base::cat("Assigning and transforming to the new CRS...\n")
        sf::st_crs(.self$LPC) <- 4326
        sf::st_crs(.self$mask) <- 4326
        sf::st_crs(.self$index) <- 4326
        .self$LPC <- sf::st_transform(.self$LPC, crs)
        .self$index <- sf::st_transform(.self$index, crs)
      }
      if (current_crs == sf::st_crs(crs)) {
        base::cat("The new CRS is the same as the current CRS. No change needed.\n")
      } else {
        base::cat("Changing and transforming to the new CRS...\n")
        sf::st_crs(.self$LPC) <- crs
        .self$LPC <- sf::st_transform(.self$LPC, crs)
        .self$index <- sf::st_transform(.self$index, crs)
      }
    },
    get_data = function() {
      return(.self$xyz)
    },
    get_lpc = function() {
      return(.self$LPC)
    },
    to_xyz = function(path) {
      utils::write.table(.self$xyz[, c("X", "Y", "Z")],
                         path,
                         row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " "
      )
    },
    to_dtm = function(resolution = 1) {
      if (is.null(.self$LPC) || !inherits(.self$LPC, "LAS")) {
        stop("LPC is not set. Initialize from a .las/.laz/.xyz before calling to_dtm().", call. = FALSE)
      }

      ground <- lidR::filter_ground(.self$LPC)
      dtm <- lidR::rasterize_terrain(ground, resolution, lidR::tin())
      .self$DTM_raw <- dtm

      if (is.null(.self$mask) || (inherits(.self$mask, "sf") && nrow(.self$mask) == 0)) {
        .self$DTM <- dtm
        return(invisible(NULL))
      }

      if (!inherits(.self$mask, "sf")) {
        stop("mask must be an sf object.", call. = FALSE)
      }
      if (any(sf::st_is_empty(sf::st_geometry(.self$mask)))) {
        stop("mask contains empty geometries.", call. = FALSE)
      }

      las_crs <- sf::st_crs(.self$LPC)
      mask_crs <- sf::st_crs(.self$mask)
      if (!is.na(las_crs) && !is.na(mask_crs) && las_crs != mask_crs) {
        .self$mask <- sf::st_transform(.self$mask, las_crs)
      }

      dtm_clip <- terra::mask(dtm, terra::vect(.self$mask))
      .self$DTM <- dtm_clip
      invisible(NULL)
    },
    to_ndsm = function(resolution = 1) {
      if (is.null(.self$LPC) || !inherits(.self$LPC, "LAS")) {
        stop("LPC is not set. Initialize from a .las/.laz/.xyz before calling to_ndsm().", call. = FALSE)
      }
      if (is.null(.self$DTM) || !inherits(.self$DTM, "SpatRaster")) {
        stop("DTM is not set. Call to_dtm() before calling to_ndsm().", call. = FALSE)
      }

      nlas <- .self$LPC - .self$DTM
      ndsmfull <- lidR::rasterize_canopy(nlas, res = resolution, lidR::p2r(0.6))
      ndsmfull[is.na(ndsmfull)] <- 0

      filled <- terra::focal(ndsmfull, w = base::matrix(1, 5, 5), fun = stats::median, na.rm = TRUE)
      clamp <- terra::clamp(filled, lower = 0)

      .self$ndsm_raw <- clamp
      .self$ndsm <- clamp
    },
    save_mask = function(path) {
      sf::st_write(.self$mask, path)
    },
    save_las = function(path) {
      lidR::writeLAS(.self$LPC, path)
    },
    save_dtm = function(path) {
      terra::writeRaster(.self$DTM, path, gdal = c("COMPRESS=LZW"), overwrite = TRUE)
    },
    save_ndsm = function(path) {
      terra::writeRaster(.self$ndsm, path, gdal = c("COMPRESS=LZW"), overwrite = TRUE)
    },
    save_sc = function(path) {
      base::save(.self, file = path)
    }
  )
)

#' Create a spatial_container
#'
#' Convenience constructor for creating a \link{spatial_container-class}
#' reference-class instance.
#'
#' @param file_path Character scalar. Path to a `.las`, `.laz`, or `.xyz` file.
#'   If empty, an empty container is created with placeholder rasters.
#'
#' @return A \link{spatial_container-class} reference-class object.
#' @export
spatial_container_new <- function(file_path = character(0)) {
  spatial_container$new(file_path = file_path)
}

