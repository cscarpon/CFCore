#' Point cloud utilities for CloudFluxCore
#'
#' Utilities for masking point clouds to a 2D polygon, denoising LAS objects
#' using SOR-based noise classification, and ICP alignment helpers.
#'
#' Notes / changes applied:
#' - Added roxygen2 blocks for each function.
#' - Namespaced non-local calls (`pkg::fun()`), including `lidR::random()` and `lidR::sor()`.
#' - Removed an unused `unique_points` object in `mask_pc()` (it had no effect).
#' - Fixed CRS handling and messaging consistency (use `message()` instead of `print()`).
#' - Avoided unqualified `LAS()` call (use `lidR::LAS()`).
#' - Avoided unqualified `projection()` call (use `lidR::projection()`).
#' - Added basic argument validation and safer temp object cleanup.
#'
#' @keywords internal

#' Create a 2D polygon mask around a point cloud
#'
#' Produces an approximate 2D footprint polygon by decimating points,
#' rasterizing their XY occupancy, polygonizing, removing holes, and simplifying.
#'
#' @param pc A `lidR::LAS` object.
#' @param decimate_keep A numeric passed to `lidR::random()` for decimation.
#'   Default keeps ~1 point per unit (see lidR docs for `random()` behavior).
#' @param raster_res Numeric raster resolution used when rasterizing occupancy.
#' @param simplify_keep Numeric in (0,1] for `rmapshaper::ms_simplify(keep=...)`.
#' @param simplify_weighting Numeric for `rmapshaper::ms_simplify(weighting=...)`.
#'
#' @return An `sf` polygon (`sf` object) representing the mask.
#' @keywords internal
mask_pc <- function(pc,
                    decimate_keep = 1,
                    raster_res = 3,
                    simplify_keep = 0.5,
                    simplify_weighting = 0.9) {
  if (is.null(pc) || !inherits(pc, "LAS")) {
    stop("`pc` must be a lidR::LAS object.", call. = FALSE)
  }

  # Step 1: reduce the number of points in the point cloud
  decimate <- lidR::decimate_points(pc, lidR::random(decimate_keep))

  # Step 2: Convert the decimated points to an sf object
  coords_sf <- sf::st_as_sf(
    decimate@data[, c("X", "Y"), drop = FALSE],
    coords = c("X", "Y"),
    crs = lidR::projection(pc)
  )

  coords_vect <- terra::vect(coords_sf)

  raster_template <- terra::rast(terra::ext(coords_vect), resolution = raster_res)

  rasterized <- terra::rasterize(coords_vect, raster_template, field = NULL, fun = "count")

  polygonized <- terra::as.polygons(rasterized, dissolve = TRUE)

  simplified_polygon <- terra::aggregate(polygonized)

  no_holes <- nngeo::st_remove_holes(sf::st_as_sf(simplified_polygon))

  final_sf <- rmapshaper::ms_simplify(
    no_holes,
    keep = simplify_keep,
    weighting = simplify_weighting,
    keep_shapes = TRUE
  )

  final_sf <- sf::st_as_sf(final_sf)

  # Ensure CRS is set and matches LAS CRS
  sf::st_crs(final_sf) <- sf::st_crs(pc)
  final_sf <- sf::st_transform(final_sf, crs = sf::st_crs(pc))

  final_sf
}

#' Denoise a LAS point cloud while removing building footprints
#'
#' Removes points within `footprint` from `area_mask` by difference, clips the LAS
#' to that resulting polygon, then runs iterative SOR noise classification and
#' removes points classified as noise (Classification == 18).
#'
#' @param laz A `lidR::LAS` object.
#' @param area_mask An `sf` polygon defining the full AOI.
#' @param footprint An `sf` polygon defining buildings to exclude.
#' @param k_sor1,m_sor1 SOR parameters for first pass (`lidR::sor(k=..., m=...)`).
#' @param k_sor2,m_sor2 SOR parameters for second pass.
#'
#' @return A `lidR::LAS` object with buildings removed and noise filtered.
#' @keywords internal
noise_filter_buildings <- function(laz, area_mask, footprint,
                                   k_sor1 = 5, m_sor1 = 3,
                                   k_sor2 = 20, m_sor2 = 5) {
  if (is.null(laz) || !inherits(laz, "LAS")) {
    stop("`laz` must be a lidR::LAS object.", call. = FALSE)
  }
  if (is.null(area_mask) || !inherits(area_mask, "sf")) {
    stop("`area_mask` must be an sf object.", call. = FALSE)
  }
  if (is.null(footprint) || !inherits(footprint, "sf")) {
    stop("`footprint` must be an sf object.", call. = FALSE)
  }

  start_time <- base::Sys.time()

  # Step 1: Create a unique PointID in the original LAS
  laz@data$PointID <- base::seq_len(base::nrow(laz@data))

  # Align CRS of inputs to LAS CRS
  laz_crs <- sf::st_crs(laz)

  if (sf::st_crs(area_mask) != laz_crs) {
    message("Transforming area_mask to LAS CRS")
    area_mask <- sf::st_transform(area_mask, crs = laz_crs)
  }

  if (sf::st_crs(footprint) != laz_crs) {
    message("Transforming footprint to LAS CRS")
    footprint <- sf::st_transform(footprint, crs = laz_crs)
  }

  diff_mask <- sf::st_difference(area_mask, footprint)

  if (sf::st_crs(diff_mask) != laz_crs) {
    message("Transforming diff_mask to LAS CRS")
    diff_mask <- sf::st_transform(diff_mask, crs = laz_crs)
  }

  message("Clipping buildings (keeping area_mask minus footprint)")
  las_no_buildings <- lidR::clip_roi(laz, diff_mask)

  message("Removing noise with iterative SOR")
  las_fix <- lidR::classify_noise(las_no_buildings, lidR::sor(k = k_sor1, m = m_sor1))
  las_fix1 <- lidR::filter_poi(las_fix, Classification != 18)
  base::rm(las_fix)

  las_fix2 <- lidR::classify_noise(las_fix1, lidR::sor(k = k_sor2, m = m_sor2))
  las_fix2 <- lidR::filter_poi(las_fix2, Classification != 18)
  base::rm(las_fix1)

  time_taken <- base::Sys.time() - start_time
  message("Denoising completed in: ", time_taken)

  las_fix2
}

#' Denoise a LAS point cloud using iterative SOR noise classification
#'
#' Runs iterative SOR noise classification and removes points classified as noise
#' (Classification == 18).
#'
#' @param laz A `lidR::LAS` object.
#' @param k_sor1,m_sor1 SOR parameters for first pass (`lidR::sor(k=..., m=...)`).
#' @param k_sor2,m_sor2 SOR parameters for second pass.
#'
#' @return A `lidR::LAS` object with noise filtered.
#' @keywords internal
noise_filter <- function(laz,
                         k_sor1 = 5, m_sor1 = 3,
                         k_sor2 = 20, m_sor2 = 5) {
  if (is.null(laz) || !inherits(laz, "LAS")) {
    stop("`laz` must be a lidR::LAS object.", call. = FALSE)
  }

  start_time <- base::Sys.time()

  laz@data$PointID <- base::seq_len(base::nrow(laz@data))

  message("Removing noise with iterative SOR")

  las_fix <- lidR::classify_noise(laz, lidR::sor(k = k_sor1, m = m_sor1))
  las_fix1 <- lidR::filter_poi(las_fix, Classification != 18)
  base::rm(las_fix)

  las_fix2 <- lidR::classify_noise(las_fix1, lidR::sor(k = k_sor2, m = m_sor2))
  las_fix2 <- lidR::filter_poi(las_fix2, Classification != 18)
  base::rm(las_fix1)

  time_taken <- base::Sys.time() - start_time
  message("Denoising completed in: ", time_taken)

  las_fix2
}

#' Transfer metadata from an original LAS to an aligned LAS by row index
#'
#' Replaces XYZ coordinates in `original_las` with those from `aligned_las`,
#' preserving all other point attributes from the original. This assumes both
#' LAS objects are in the same point order. If voxelization or reindexing was
#' applied, this will not be a valid mapping.
#'
#' @param original_las A `lidR::LAS` object to receive updated XYZ.
#' @param aligned_las A `lidR::LAS` object providing updated XYZ.
#' @param crs CRS to assign to the returned LAS (e.g., `sf::st_crs(original_las)`).
#'
#' @return A `lidR::LAS` with original attributes and aligned XYZ.
#' @keywords internal
pc_metadata <- function(original_las, aligned_las, crs) {
  if (is.null(original_las) || is.null(aligned_las)) {
    stop("Original or aligned LAS object is NULL.", call. = FALSE)
  }
  if (!inherits(original_las, "LAS") || !inherits(aligned_las, "LAS")) {
    stop("`original_las` and `aligned_las` must be lidR::LAS objects.", call. = FALSE)
  }

  original_data <- original_las@data
  aligned_data <- aligned_las@data

  min_len <- base::min(base::nrow(original_data), base::nrow(aligned_data))
  original_data <- original_data[1:min_len, , drop = FALSE]
  aligned_data <- aligned_data[1:min_len, , drop = FALSE]

  original_data$X <- aligned_data$X
  original_data$Y <- aligned_data$Y
  original_data$Z <- aligned_data$Z

  original_las@data <- original_data

  sf::st_crs(original_las) <- crs

  original_las
}

#' Align two LAS point clouds using ICP after voxelization
#'
#' Voxelizes both point clouds, runs `Morpho::icpmat()` to align target to source,
#' and writes an aligned LAS to disk. Metadata is transferred from the voxelized
#' target points via `point_source_id`.
#'
#' Important: This returns a voxel-level LAS, not a full-resolution aligned LAS.
#' If you need full resolution alignment, you should apply the resulting
#' transformation to the original target point matrix.
#'
#' @param source_las `lidR::LAS` reference point cloud (fixed).
#' @param target_las `lidR::LAS` point cloud to align to the source (moving).
#' @param output_las_path Output path for `lidR::writeLAS()`.
#' @param voxel_size Numeric voxel size passed to `lidR::voxelize_points()`.
#' @param max_iter Integer ICP iterations passed to `Morpho::icpmat()`.
#' @param threads Integer thread count passed to `Morpho::icpmat()`.
#'
#' @return A voxelized `lidR::LAS` object representing the aligned target.
#' @keywords internal
align_las_icp_voxelized <- function(source_las, target_las, output_las_path,
                                    voxel_size = 0.2,
                                    max_iter = 50,
                                    threads = 12) {
  if (is.null(source_las) || !inherits(source_las, "LAS")) {
    stop("`source_las` must be a lidR::LAS object.", call. = FALSE)
  }
  if (is.null(target_las) || !inherits(target_las, "LAS")) {
    stop("`target_las` must be a lidR::LAS object.", call. = FALSE)
  }
  if (!is.character(output_las_path) || base::length(output_las_path) != 1L) {
    stop("`output_las_path` must be a single character path.", call. = FALSE)
  }

  source_las@data$PointID <- base::seq_len(base::nrow(source_las@data))
  target_las@data$PointID <- base::seq_len(base::nrow(target_las@data))

  target_metadata <- target_las@data

  las_source_voxel <- lidR::voxelize_points(source_las, voxel_size)
  las_target_voxel <- lidR::voxelize_points(target_las, voxel_size)

  if (is.null(las_target_voxel@data$point_source_id)) {
    stop("Voxelized target LAS does not contain `point_source_id`. Cannot map metadata.", call. = FALSE)
  }

  idx <- las_target_voxel@data$point_source_id
  if (base::any(idx < 1L) || base::any(idx > base::nrow(target_metadata))) {
    stop("`point_source_id` indices are out of bounds for target metadata.", call. = FALSE)
  }

  voxel_metadata <- target_metadata[idx, , drop = FALSE]

  source_pts <- base::as.matrix(las_source_voxel@data[, c("X", "Y", "Z"), drop = FALSE])
  target_pts <- base::as.matrix(las_target_voxel@data[, c("X", "Y", "Z"), drop = FALSE])

  icp_result <- Morpho::icpmat(target_pts, source_pts, iterations = max_iter, threads = threads)

  aligned_pts <- (icp_result$R %*% base::t(target_pts)) + icp_result$t
  aligned_pts <- base::t(aligned_pts)

  aligned_las_data <- base::data.frame(
    X = aligned_pts[, 1],
    Y = aligned_pts[, 2],
    Z = aligned_pts[, 3]
  )

  drop_cols <- base::names(voxel_metadata) %in% c("X", "Y", "Z")
  aligned_las_data <- base::cbind(aligned_las_data, voxel_metadata[, !drop_cols, drop = FALSE])

  aligned_las <- lidR::LAS(aligned_las_data)

  lidR::writeLAS(aligned_las, output_las_path)

  message("ICP alignment with voxelization complete. Output saved to: ", output_las_path)

  aligned_las
}
