#' cloudFlux reference class
#'
#' A mutable workflow container that orchestrates two-epoch point cloud change
#' processing in CFCore: mask generation, denoising, DTM/nDSM derivation, raster
#' alignment, and continuous difference computation.
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
    icp_aligned_path = "character",
    use_icp = "logical",
    use_gpu = "logical",
    voxel_size = "numeric",
    icp_method = "character"
  ),
  methods = list(

    initialize = function(source_path = character(0),
                          target_path = character(0),
                          source_las = NULL,
                          target_las = NULL,
                          epsg = 26917,
                          resolution = 1,
                          use_icp = FALSE,
                          use_gpu = TRUE,
                          voxel_size = 0.05,
                          icp_method = "point-to-plane") {

      .self$epsg <- base::as.integer(epsg)
      .self$resolution <- resolution
      .self$timings <- list()
      .self$icp_aligned_path <- character(0)

      .self$use_icp <- use_icp
      .self$use_gpu <- use_gpu
      .self$voxel_size <- voxel_size
      .self$icp_method <- icp_method

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

      invisible(.self)
    },

    run = function(resolution = .self$resolution,
                   method = "bilinear",
                   k_sor1 = 5, m_sor1 = 3,
                   k_sor2 = 20, m_sor2 = 5) {

      message("Setting coordinate reference systems...")
      if (!is.null(.self$pc_source$LPC)) .self$pc_source$set_crs(.self$epsg)
      if (!is.null(.self$pc_target$LPC)) .self$pc_target$set_crs(.self$epsg)

      if (isTRUE(.self$use_icp)) {
        message(sprintf("Running ICP Alignment (GPU requested: %s)...", .self$use_gpu))
        .self$icp_align_open3d(
          voxel_size = .self$voxel_size,
          icp_method = .self$icp_method,
          use_gpu = .self$use_gpu
        )
      }

      message("Building masks...")
      .self$build_masks()

      message("Applying noise filters...")
      .self$denoise(k_sor1 = k_sor1, m_sor1 = m_sor1, k_sor2 = k_sor2, m_sor2 = m_sor2)

      message("Generating DTM and nDSM rasters...")
      .self$build_rasters(resolution = resolution)

      message("Aligning rasters...")
      .self$align_rasters(method = method)

      message("Computing continuous differences...")
      .self$compute_diffs()

      message("Workflow complete.")
      invisible(.self)
    },

    run_all = function(...) {
      .self$run(...)
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
        .self$pc_source$LPC, k_sor1 = k_sor1, m_sor1 = m_sor1, k_sor2 = k_sor2, m_sor2 = m_sor2
      )
      .self$pc_target$LPC <- noise_filter(
        .self$pc_target$LPC, k_sor1 = k_sor1, m_sor1 = m_sor1, k_sor2 = k_sor2, m_sor2 = m_sor2
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

    icp_align_open3d = function(icp_py = "inst/py/icp_hybrid.py", # Make sure this matches your script name!
                                icp_condaenv = "icp_conda",
                                voxel_size = 0.05,
                                icp_method = "point-to-plane",
                                use_gpu = TRUE) {
      if (!base::requireNamespace("reticulate", quietly = TRUE)) stop("Install 'reticulate'.")
      if (!base::requireNamespace("lidR", quietly = TRUE)) stop("Install 'lidR'.")
      if (is.null(.self$pc_source$LPC) || is.null(.self$pc_target$LPC)) stop("Clouds must be set.")

      if (!base::file.exists(icp_py)) {
        icp_py2 <- base::system.file("py", "icp_hybrid.py", package = "CFCore")
        if (base::nzchar(icp_py2)) icp_py <- icp_py2
      }
      if (!base::file.exists(icp_py)) stop("ICP python script not found: ", icp_py)

      t0 <- base::Sys.time() # START TIMER

      source_path <- base::file.path(base::tempdir(), "cfcore_icp_source.laz")
      target_path <- base::file.path(base::tempdir(), "cfcore_icp_target.laz")

      lidR::writeLAS(.self$pc_source$LPC, source_path)
      lidR::writeLAS(.self$pc_target$LPC, target_path)

      reticulate::use_condaenv(icp_condaenv, required = TRUE)
      reticulate::source_python(icp_py)

      icp_aligner <- HybridICP(
        source_path, target_path,
        voxel_size = voxel_size,
        icp_method = icp_method,
        use_gpu = use_gpu
      )

      res <- icp_aligner$align()

      aligned_path <- tryCatch(res[[1]], error = function(e) NULL)
      msg <- tryCatch(res[[2]], error = function(e) NULL)

      if (is.null(aligned_path) || !base::nzchar(base::as.character(aligned_path))) {
        stop("ICP failed to return path. Msg: ", base::as.character(msg), call. = FALSE)
      }

      .self$icp_aligned_path <- base::as.character(aligned_path)
      .self$timings$icp_message <- base::as.character(msg)

      .self$pc_target <- spatial_container$new(file_path = .self$icp_aligned_path)
      .self$pc_target$set_crs(.self$epsg)

      # STOP TIMER AND PRINT
      t_end <- base::Sys.time()
      duration <- base::round(base::difftime(t_end, t0, units = "secs"), 2)
      message(sprintf("ICP Alignment completed in %s seconds.", duration))

      .self$timings$icp_total <- duration
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
        stop("Rasters must be aligned before computing diffs.", call. = FALSE)
      }

      src_ndsm_fixed <- terra::subst(.self$aligned_ndsm$source, NA, 0)
      tgt_ndsm_fixed <- terra::subst(.self$aligned_ndsm$target, NA, 0)

      .self$ndsm_diff <- diff_ndsm(src_ndsm_fixed, tgt_ndsm_fixed)
      .self$dtm_diff  <- diff_dtm(.self$aligned_dtm$source, .self$aligned_dtm$target)

      if (!is.null(.self$mask_union)) {
        .self$ndsm_diff <- terra::mask(.self$ndsm_diff, .self$mask_union)
        .self$dtm_diff  <- terra::mask(.self$dtm_diff,  .self$mask_union)
      }

      .self$timings$diffs <- base::Sys.time() - t0
      invisible(.self)
    },

    summary_stats = function() {
      list(
        ndsm = diff_summary_stats(.self$ndsm_diff),
        dtm  = diff_summary_stats(.self$dtm_diff)
      )
    },

    summary_by_breaks = function(breaks, unit = "m") {
      list(
        ndsm = diff_summary_by_breaks(.self$ndsm_diff, breaks, unit = unit),
        dtm  = diff_summary_by_breaks(.self$dtm_diff,  breaks, unit = unit)
      )
    },

    plot_hist = function(which = "ndsm", bins = 60, title = NULL) {
      which <- base::tolower(which)
      r <- if (identical(which, "dtm")) .self$dtm_diff else .self$ndsm_diff
      if (is.null(title)) title <- if (identical(which, "dtm")) "DTM difference" else "nDSM difference"
      plot_diff_hist(r, bins = bins, title = title)
    },

    plot_binned = function(breaks, which = "ndsm", title = NULL) {
      which <- base::tolower(which)
      r <- if (identical(which, "dtm")) .self$dtm_diff else .self$ndsm_diff
      if (is.null(title)) title <- if (identical(which, "dtm")) "DTM difference" else "nDSM difference"
      plot_diff_binned(r, breaks = breaks, title = title)
    },

    # --- UPDATED QUANTILE DEFAULT (c(0, 1)) TO PREVENT HOLES IN LEAFLET MAPS ---
    map = function(diff_breaks = NULL, epsg = 4326, diff_palette = "viridis",
                   diff_range = NULL, diff_quantiles = c(0, 1), diff_break_colors = NULL) {
      displayMap(
        dtm1 = .self$pc_source$DTM, ndsm1 = .self$pc_source$ndsm,
        dtm2 = .self$pc_target$DTM, ndsm2 = .self$pc_target$ndsm,
        dtm_diff = .self$dtm_diff, ndsm_diff = .self$ndsm_diff,
        area_mask = .self$mask_union, epsg = epsg, diff_palette = diff_palette,
        diff_range = diff_range, diff_quantiles = diff_quantiles,
        diff_breaks = diff_breaks, diff_break_colors = diff_break_colors
      )
    },

    save_rasters = function(out_dir = ".", prefix = "cf") {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

      if (!is.null(.self$ndsm_diff)) {
        terra::writeRaster(.self$ndsm_diff, file.path(out_dir, paste0(prefix, "_ndsm_diff.tif")), overwrite = TRUE)
      }
      if (!is.null(.self$dtm_diff)) {
        terra::writeRaster(.self$dtm_diff, file.path(out_dir, paste0(prefix, "_dtm_diff.tif")), overwrite = TRUE)
      }
      if (!is.null(.self$aligned_dtm$source)) {
        terra::writeRaster(.self$aligned_dtm$source, file.path(out_dir, paste0(prefix, "_dtm_source.tif")), overwrite = TRUE)
        terra::writeRaster(.self$aligned_dtm$target, file.path(out_dir, paste0(prefix, "_dtm_target.tif")), overwrite = TRUE)
      }
      if (!is.null(.self$aligned_ndsm$source)) {
        terra::writeRaster(.self$aligned_ndsm$source, file.path(out_dir, paste0(prefix, "_ndsm_source.tif")), overwrite = TRUE)
        terra::writeRaster(.self$aligned_ndsm$target, file.path(out_dir, paste0(prefix, "_ndsm_target.tif")), overwrite = TRUE)
      }
      message("Rasters successfully saved to ", base::normalizePath(out_dir))
      invisible(.self)
    }
  )
)

#' Create a cloudFlux object
#'
#' Convenience constructor for the cloudFlux reference class workflow.
#'
#' @param source_path Character path to the source LAS/LAZ file.
#' @param target_path Character path to the target LAS/LAZ file.
#' @param source_las `lidR::LAS` object for the source epoch.
#' @param target_las `lidR::LAS` object for the target epoch.
#' @param epsg Integer EPSG code for the working projection (default: 26917).
#' @param resolution Numeric raster resolution for DTM/nDSM (default: 1).
#' @param use_icp Logical; if TRUE, runs ICP alignment automatically (default: FALSE).
#' @param use_gpu Logical; if TRUE, attempts GPU acceleration for ICP (default: TRUE).
#' @param voxel_size Numeric voxel size for ICP downsampling (default: 0.05).
#' @param icp_method Character; ICP method, e.g., "point-to-plane" (default).
#'
#' @import reticulate
#' @return A `cloudFlux` reference class object.
#' @export
cloudFlux_new <- function(source_path = character(0),
                          target_path = character(0),
                          source_las = NULL,
                          target_las = NULL,
                          epsg = 26917,
                          resolution = 1,
                          use_icp = FALSE,
                          use_gpu = TRUE,
                          voxel_size = 0.05,
                          icp_method = "point-to-plane") {
  cf <- cloudFlux$new(
    source_path = source_path,
    target_path = target_path,
    source_las = source_las,
    target_las = target_las,
    epsg = epsg,
    resolution = resolution,
    use_icp = use_icp,
    use_gpu = use_gpu,
    voxel_size = voxel_size,
    icp_method = icp_method
  )
  return(cf)
}
