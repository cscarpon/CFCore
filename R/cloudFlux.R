#' cloudFlux reference class
#'
#' A mutable workflow container that orchestrates two-epoch point cloud change
#' processing in CFCore: mask generation, denoising, DTM/nDSM derivation, raster
#' alignment, and continuous difference computation.
#'
#' This class is designed to reduce boilerplate for common two-epoch workflows.
#' For API-style usage, you can still call the underlying standalone functions
#' directly (`mask_pc()`, `noise_filter()`, `process_raster()`, `diff_ndsm()`, etc.).
#'
#' @section Fields:
#' \describe{
#'   \item{epsg}{`integer`. Working EPSG used for point clouds and masks.}
#'   \item{resolution}{`numeric`. Raster resolution used for DTM/nDSM.}
#'   \item{pc_source}{`spatial_container`. Source (earlier) epoch container.}
#'   \item{pc_target}{`spatial_container`. Target (later) epoch container.}
#'   \item{mask_union}{`terra::SpatVector`. Union mask used for raster masking/alignment.}
#'   \item{aligned_ndsm}{`list`. Output of `process_raster()` for nDSM rasters.}
#'   \item{aligned_dtm}{`list`. Output of `process_raster()` for DTM rasters.}
#'   \item{ndsm_diff}{`terra::SpatRaster`. Continuous nDSM difference (target - source).}
#'   \item{dtm_diff}{`terra::SpatRaster`. Continuous DTM difference (target - source).}
#'   \item{timings}{`list`. Named timings captured during processing (optional).}
#'   \item{icp_aligned_path}{`character`. Path to ICP-aligned output (if used).}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{initialize(...)}{Optionally loads LAS/LAZ, sets CRS, and can run the workflow.}
#'   \item{run(...)}{Runs: masks, denoise, DTM, nDSM, align rasters, compute diffs.}
#'   \item{run_all(...)}{Alias for `run()` (kept for compatibility).}
#'   \item{build_masks()}{Computes 2D masks using `mask_pc()`.}
#'   \item{denoise()}{Runs `noise_filter()` on both epochs.}
#'   \item{build_rasters(resolution = self$resolution)}{Calls `to_dtm()` and `to_ndsm()` on both epochs.}
#'   \item{icp_align_open3d(...)}{Optional Open3D ICP alignment via reticulate (aligns target to source).}
#'   \item{align_rasters(method = "bilinear")}{Aligns DTM and nDSM rasters using `process_raster()`.}
#'   \item{compute_diffs()}{Computes continuous differences and masks them to the union mask.}
#'   \item{summary_stats()}{Returns `list(ndsm=..., dtm=...)` from `diff_summary_stats()`.}
#'   \item{summary_by_breaks(breaks, unit = "m")}{Returns binned summaries via `diff_summary_by_breaks()`.}
#'   \item{plot_hist(which = "ndsm", bins = 60, title = NULL)}{Histogram plot for a diff raster.}
#'   \item{plot_binned(breaks, which = "ndsm", title = NULL)}{Binned bar plot for a diff raster.}
#'   \item{map(...)}{Leaflet map via `displayMap()` using stored rasters/diffs.}
#' }
#'
#' @name cloudFlux
#' @docType class
#' @aliases cloudFlux-class
#' @keywords internal
NULL

# Avoid hard stops during package load; enforce at runtime inside methods instead.
methods::setOldClass("sf")
methods::setOldClass("sfc")

#' @import methods
#' @importClassesFrom lidR LAS
#' @importClassesFrom terra SpatRaster SpatVector
cloudFlux <- methods::setRefClass(
  "cloudFlux",
  fields = list(
    epsg = "numeric",
    resolution = "numeric",
    pc_source = "ANY",
    pc_target = "ANY",
    mask_union = "ANY",
    aligned_ndsm = "list",
    aligned_dtm = "list",
    ndsm_diff = "ANY",
    dtm_diff = "ANY",
    timings = "list",
    icp_aligned_path = "character"
  ),
  methods = list(

    initialize = function(source_path = character(0),
                          target_path = character(0),
                          source_las = NULL,
                          target_las = NULL,
                          epsg = 26917,
                          resolution = 1,
                          run = FALSE,
                          icp_align = FALSE,
                          icp_py = "inst/py/icp_open3D.py",
                          icp_condaenv = "icp_conda",
                          voxel_size = 0.05,
                          icp_method = "point-to-plane",
                          method = "bilinear",
                          k_sor1 = 5, m_sor1 = 3,
                          k_sor2 = 20, m_sor2 = 5) {

      .self$epsg <- base::as.integer(epsg)
      .self$resolution <- resolution
      .self$timings <- list()
      .self$icp_aligned_path <- character(0)

      .self$pc_source <- spatial_container$new()
      .self$pc_target <- spatial_container$new()

      if (!is.null(source_las)) {
        .self$pc_source$LPC <- source_las
      } else if (base::length(source_path) > 0) {
        .self$pc_source <- spatial_container$new(file_path = source_path)
      }

      if (!is.null(target_las)) {
        .self$pc_target$LPC <- target_las
      } else if (base::length(target_path) > 0) {
        .self$pc_target <- spatial_container$new(file_path = target_path)
      }

      if (!is.null(.self$pc_source$LPC)) .self$pc_source$set_crs(.self$epsg)
      if (!is.null(.self$pc_target$LPC)) .self$pc_target$set_crs(.self$epsg)

      if (isTRUE(icp_align)) {
        .self$icp_align_open3d(
          icp_py = icp_py,
          icp_condaenv = icp_condaenv,
          voxel_size = voxel_size,
          icp_method = icp_method
        )
      }

      if (isTRUE(run)) {
        .self$run(
          resolution = resolution,
          method = method,
          k_sor1 = k_sor1, m_sor1 = m_sor1,
          k_sor2 = k_sor2, m_sor2 = m_sor2
        )
      }

      invisible(.self)
    },

    run = function(resolution = .self$resolution,
                   method = "bilinear",
                   k_sor1 = 5, m_sor1 = 3,
                   k_sor2 = 20, m_sor2 = 5) {
      .self$build_masks()
      .self$denoise(k_sor1 = k_sor1, m_sor1 = m_sor1, k_sor2 = k_sor2, m_sor2 = m_sor2)
      .self$build_rasters(resolution = resolution)
      .self$align_rasters(method = method)
      .self$compute_diffs()
      invisible(.self)
    },

    # compatibility alias
    run_all = function(resolution = .self$resolution,
                       method = "bilinear",
                       k_sor1 = 5, m_sor1 = 3,
                       k_sor2 = 20, m_sor2 = 5) {
      .self$run(
        resolution = resolution, method = method,
        k_sor1 = k_sor1, m_sor1 = m_sor1,
        k_sor2 = k_sor2, m_sor2 = m_sor2
      )
      invisible(.self)
    },

    build_masks = function() {
      t0 <- base::Sys.time()
      if (is.null(.self$pc_source$LPC) || is.null(.self$pc_target$LPC)) {
        stop("Both pc_source$LPC and pc_target$LPC must be set before building masks.", call. = FALSE)
      }
      .self$pc_source$mask <- mask_pc(.self$pc_source$LPC)
      .self$pc_target$mask <- mask_pc(.self$pc_target$LPC)
      .self$timings$mask <- base::Sys.time() - t0
      invisible(.self)
    },

    denoise = function(k_sor1 = 5, m_sor1 = 3, k_sor2 = 20, m_sor2 = 5) {
      t0 <- base::Sys.time()
      .self$pc_source$LPC <- noise_filter(
        .self$pc_source$LPC,
        k_sor1 = k_sor1, m_sor1 = m_sor1,
        k_sor2 = k_sor2, m_sor2 = m_sor2
      )
      .self$pc_target$LPC <- noise_filter(
        .self$pc_target$LPC,
        k_sor1 = k_sor1, m_sor1 = m_sor1,
        k_sor2 = k_sor2, m_sor2 = m_sor2
      )
      .self$timings$denoise <- base::Sys.time() - t0
      invisible(.self)
    },

    build_rasters = function(resolution = .self$resolution) {
      t0 <- base::Sys.time()
      .self$pc_source$to_dtm(resolution = resolution)
      .self$pc_target$to_dtm(resolution = resolution)
      .self$pc_source$to_ndsm(resolution = resolution)
      .self$pc_target$to_ndsm(resolution = resolution)
      .self$timings$rasters <- base::Sys.time() - t0
      invisible(.self)
    },

    icp_align_open3d = function(icp_py = "inst/py/icp_open3D.py",
                                icp_condaenv = "icp_conda",
                                voxel_size = 0.05,
                                icp_method = "point-to-plane") {
      if (!base::requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required for ICP alignment but is not installed.", call. = FALSE)
      }
      if (!base::requireNamespace("lidR", quietly = TRUE)) {
        stop("Package 'lidR' is required for ICP alignment but is not installed.", call. = FALSE)
      }

      if (is.null(.self$pc_source$LPC) || is.null(.self$pc_target$LPC)) {
        stop("Both pc_source$LPC and pc_target$LPC must be set before ICP alignment.", call. = FALSE)
      }

      # Resolve icp_py robustly for installed package or dev workflows
      if (!base::file.exists(icp_py)) {
        icp_py2 <- base::system.file("py", "icp_open3D.py", package = "CFCore")
        if (base::nzchar(icp_py2)) icp_py <- icp_py2
      }
      if (!base::file.exists(icp_py)) {
        stop("ICP python script not found at: ", icp_py, call. = FALSE)
      }

      t0 <- base::Sys.time()

      source_path <- base::file.path(base::tempdir(), "cfcore_icp_source.laz")
      target_path <- base::file.path(base::tempdir(), "cfcore_icp_target.laz")

      lidR::writeLAS(.self$pc_source$LPC, source_path)
      lidR::writeLAS(.self$pc_target$LPC, target_path)

      reticulate::use_condaenv(icp_condaenv, required = TRUE)
      reticulate::source_python(icp_py)

      if (!base::exists("Open3DICP")) {
        stop("Open3DICP was not loaded from the python script. Check icp_py exports.", call. = FALSE)
      }

      icp_aligner <- Open3DICP(
        source_path, target_path,
        voxel_size = voxel_size,
        icp_method = icp_method
      )

      res <- icp_aligner$align()

      # Python align() returns (aligned_path, message). Reticulate maps to list/tuple.
      aligned_path <- NULL
      msg <- NULL

      if (base::is.list(res) && base::length(res) >= 1L) {
        aligned_path <- res[[1]]
        if (base::length(res) >= 2L) msg <- res[[2]]
      } else {
        aligned_path <- tryCatch(res[[1]], error = function(e) NULL)
        msg <- tryCatch(res[[2]], error = function(e) NULL)
      }

      if (is.null(aligned_path) || !base::nzchar(base::as.character(aligned_path))) {
        stop("ICP did not return an aligned path. Message: ", base::as.character(msg), call. = FALSE)
      }

      aligned_path <- base::as.character(aligned_path)

      if (!base::file.exists(aligned_path)) {
        stop(
          "ICP returned a path that does not exist: ", aligned_path,
          "\nMessage: ", base::as.character(msg),
          call. = FALSE
        )
      }

      .self$icp_aligned_path <- aligned_path
      .self$timings$icp_message <- base::as.character(msg)

      # Replace target point cloud with aligned version
      .self$pc_target <- spatial_container$new(file_path = .self$icp_aligned_path)
      .self$pc_target$set_crs(.self$epsg)

      .self$timings$icp_open3d <- base::Sys.time() - t0
      invisible(.self)
    },

    align_rasters = function(method = "bilinear") {
      t0 <- base::Sys.time()
      if (is.null(.self$pc_source$mask) || is.null(.self$pc_target$mask)) {
        stop("Masks must be created before raster alignment. Run build_masks() first.", call. = FALSE)
      }

      .self$aligned_ndsm <- process_raster(
        source = .self$pc_source$ndsm_raw,
        target = .self$pc_target$ndsm_raw,
        source_mask = .self$pc_source$mask,
        target_mask = .self$pc_target$mask,
        method = method
      )

      .self$aligned_dtm <- process_raster(
        source = .self$pc_source$DTM_raw,
        target = .self$pc_target$DTM_raw,
        source_mask = .self$pc_source$mask,
        target_mask = .self$pc_target$mask,
        method = method
      )

      .self$mask_union <- .self$aligned_ndsm$vect_mask
      .self$timings$align_rasters <- base::Sys.time() - t0
      invisible(.self)
    },

    compute_diffs = function() {
      t0 <- base::Sys.time()
      if (base::length(.self$aligned_ndsm) == 0L || base::length(.self$aligned_dtm) == 0L) {
        stop("Rasters must be aligned before computing diffs. Run align_rasters() first.", call. = FALSE)
      }

      .self$ndsm_diff <- diff_ndsm(.self$aligned_ndsm$source, .self$aligned_ndsm$target)
      .self$dtm_diff  <- diff_dtm(.self$aligned_dtm$source,  .self$aligned_dtm$target)

      if (!is.null(.self$mask_union)) {
        .self$ndsm_diff <- terra::mask(.self$ndsm_diff, .self$mask_union)
        .self$dtm_diff  <- terra::mask(.self$dtm_diff,  .self$mask_union)
      }

      .self$timings$diffs <- base::Sys.time() - t0
      invisible(.self)
    },

    summary_stats = function() {
      if (is.null(.self$ndsm_diff) || is.null(.self$dtm_diff)) {
        stop("Diff rasters are not available. Run `run()` first.", call. = FALSE)
      }
      list(
        ndsm = diff_summary_stats(.self$ndsm_diff),
        dtm  = diff_summary_stats(.self$dtm_diff)
      )
    },

    summary_by_breaks = function(breaks, unit = "m") {
      if (is.null(.self$ndsm_diff) || is.null(.self$dtm_diff)) {
        stop("Diff rasters are not available. Run `run()` first.", call. = FALSE)
      }
      list(
        ndsm = diff_summary_by_breaks(.self$ndsm_diff, breaks, unit = unit),
        dtm  = diff_summary_by_breaks(.self$dtm_diff,  breaks, unit = unit)
      )
    },

    plot_hist = function(which = "ndsm", bins = 60, title = NULL) {
      if (is.null(.self$ndsm_diff) || is.null(.self$dtm_diff)) {
        stop("Diff rasters are not available. Run `run()` first.", call. = FALSE)
      }
      which <- base::tolower(which)
      r <- if (identical(which, "dtm")) .self$dtm_diff else .self$ndsm_diff
      if (is.null(title)) title <- if (identical(which, "dtm")) "DTM difference histogram" else "nDSM difference histogram"
      plot_diff_hist(r, bins = bins, title = title)
    },

    plot_binned = function(breaks, which = "ndsm", title = NULL) {
      if (is.null(.self$ndsm_diff) || is.null(.self$dtm_diff)) {
        stop("Diff rasters are not available. Run `run()` first.", call. = FALSE)
      }
      which <- base::tolower(which)
      r <- if (identical(which, "dtm")) .self$dtm_diff else .self$ndsm_diff
      if (is.null(title)) title <- if (identical(which, "dtm")) "DTM difference by bins" else "nDSM difference by bins"
      plot_diff_binned(r, breaks = breaks, title = title)
    },

    map = function(diff_breaks = NULL,
                   epsg = 4326,
                   diff_palette = "viridis",
                   diff_range = NULL,
                   diff_quantiles = c(0.02, 0.98),
                   diff_break_colors = NULL) {
      if (is.null(.self$pc_source$DTM) || is.null(.self$pc_target$DTM) ||
          is.null(.self$pc_source$ndsm) || is.null(.self$pc_target$ndsm)) {
        stop("Base rasters are not available. Run `run()` first.", call. = FALSE)
      }

      displayMap(
        dtm1 = .self$pc_source$DTM, ndsm1 = .self$pc_source$ndsm,
        dtm2 = .self$pc_target$DTM, ndsm2 = .self$pc_target$ndsm,
        dtm_diff = .self$dtm_diff, ndsm_diff = .self$ndsm_diff,
        area_mask = .self$mask_union,
        epsg = epsg,
        diff_palette = diff_palette,
        diff_range = diff_range,
        diff_quantiles = diff_quantiles,
        diff_breaks = diff_breaks,
        diff_break_colors = diff_break_colors
      )
    }
  )
)

#' Create a cloudFlux object
#'
#' Convenience constructor for the cloudFlux reference class.
#'
#' @param source_path,target_path Character paths to LAS/LAZ inputs (optional if LAS provided).
#' @param source_las,target_las `lidR::LAS` objects (optional if paths provided).
#' @param epsg Integer EPSG used for the workflow (e.g., 26917).
#' @param resolution Numeric raster resolution used for DTM/nDSM.
#' @param run Logical; if TRUE runs the full workflow after initialization.
#' @param icp_align Logical; if TRUE runs Open3D ICP (aligns target to source) before raster steps.
#' @param icp_py Path to the Open3D ICP python script.
#' @param icp_condaenv Conda environment name for reticulate to use.
#' @param voxel_size,icp_method ICP parameters passed into the python class.
#' @param method Resampling method for `process_raster()` (default "bilinear").
#' @param k_sor1,m_sor1,k_sor2,m_sor2 SOR noise filter parameters.
#'
#' @return A cloudFlux reference-class object.
#' @export
cloudFlux_new <- function(source_path = character(0),
                          target_path = character(0),
                          source_las = NULL,
                          target_las = NULL,
                          epsg = 26917,
                          resolution = 1,
                          run = FALSE,
                          icp_align = FALSE,
                          icp_py = "inst/py/icp_open3D.py",
                          icp_condaenv = "icp_conda",
                          voxel_size = 0.05,
                          icp_method = "point-to-plane",
                          method = "bilinear",
                          k_sor1 = 5, m_sor1 = 3,
                          k_sor2 = 20, m_sor2 = 5) {
  cloudFlux$new(
    source_path = source_path,
    target_path = target_path,
    source_las = source_las,
    target_las = target_las,
    epsg = epsg,
    resolution = resolution,
    run = run,
    icp_align = icp_align,
    icp_py = icp_py,
    icp_condaenv = icp_condaenv,
    voxel_size = voxel_size,
    icp_method = icp_method,
    method = method,
    k_sor1 = k_sor1, m_sor1 = m_sor1,
    k_sor2 = k_sor2, m_sor2 = m_sor2
  )
}
