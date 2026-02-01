#' Raster and polygon utilities for CloudFluxCore
#'
#' Internal utilities for CRS alignment, masking, resampling, continuous
#' difference computation, and optional reporting summaries/plots.
#'
#' Design:
#' - Core outputs are continuous rasters (API-friendly).
#' - Binning (breaks) is optional and used only for reporting/plots/legends.
#'
#' @keywords internal

#' Check whether an sfc geometry object is empty
#'
#' @param sfc An `sf::sfc` object.
#' @return `TRUE` if length is zero, otherwise `FALSE`.
#' @keywords internal
is_empty_sfc <- function(sfc) {
  base::length(sfc) == 0L
}

#' Transform a polygon CRS to match a target polygon CRS
#'
#' If `source_polygon` has missing CRS, assigns the provided `crs` before
#' transforming to match `target_polygon`.
#'
#' @param source_polygon `sf` object.
#' @param target_polygon `sf` object whose CRS is treated as authoritative.
#' @param crs EPSG code (integer or coercible).
#' @return `source_polygon` transformed to the CRS of `target_polygon`.
#' @keywords internal
transform_polygon_crs <- function(source_polygon, target_polygon, crs) {
  crs_int <- base::as.integer(crs)

  if (is.na(sf::st_crs(source_polygon))) {
    sf::st_crs(source_polygon) <- sf::st_crs(crs_int)
  }

  if (sf::st_crs(source_polygon) != sf::st_crs(target_polygon)) {
    source_polygon <- sf::st_transform(source_polygon, sf::st_crs(target_polygon))
  }

  source_polygon
}

#' Transform a raster CRS to match a target raster CRS
#'
#' If `source_raster` has missing CRS, assigns the provided `crs` (EPSG) before
#' projecting to match `target_raster`.
#'
#' @param source_raster `terra::SpatRaster`.
#' @param target_raster `terra::SpatRaster` whose CRS is treated as authoritative.
#' @param crs EPSG code (integer or coercible) used only when `source_raster` CRS is missing.
#' @return `source_raster` projected to match `target_raster` if CRS differs.
#' @keywords internal
transform_raster_crs <- function(source_raster, target_raster, crs) {
  crs_input <- base::paste0("EPSG:", base::as.integer(crs))

  if (is.na(terra::crs(source_raster))) {
    terra::crs(source_raster) <- crs_input
  }

  if (terra::crs(source_raster) != terra::crs(target_raster)) {
    source_raster <- terra::project(source_raster, target_raster)
  }

  source_raster
}

#' Frequency table with percentages for a classified raster
#'
#' This is only meaningful when `raster` is categorical/classified. For continuous
#' rasters, use `diff_summary_stats()` or `diff_summary_by_breaks()`.
#'
#' @param raster `terra::SpatRaster` categorical raster.
#' @return A data.frame with counts and percentages by class value. If the raster
#'   has no frequencies (e.g., all NA), returns an empty data.frame with columns
#'   `count` and `percentage`.
#' @keywords internal
diff_values <- function(raster) {
  class_freq <- terra::freq(raster)

  if (is.null(class_freq) || base::nrow(class_freq) == 0L) {
    return(base::data.frame(count = integer(0), percentage = numeric(0)))
  }

  total_cells <- base::sum(class_freq[, "count"])
  class_freq$percentage <- (class_freq$count / total_cells) * 100

  base::rownames(class_freq) <- class_freq$value
  class_freq <- class_freq[, -1, drop = FALSE]

  class_freq
}

#' Align, resample, and mask a raster pair to a shared union mask
#'
#' @param source `terra::SpatRaster`.
#' @param target `terra::SpatRaster`.
#' @param source_mask `sf` polygon(s) representing the source AOI.
#' @param target_mask `sf` polygon(s) representing the target AOI.
#' @param method Character. Resampling method passed to `terra::resample()`.
#' @return A list with `source`, `target`, and `vect_mask` (a `terra::SpatVector`).
#' @keywords internal
process_raster <- function(source, target, source_mask, target_mask, method = "bilinear") {
  if (!inherits(source, "SpatRaster") || !inherits(target, "SpatRaster")) {
    stop("`source` and `target` must be terra::SpatRaster.", call. = FALSE)
  }
  if (!inherits(source_mask, "sf") || !inherits(target_mask, "sf")) {
    stop("`source_mask` and `target_mask` must be sf objects.", call. = FALSE)
  }

  # Ensure same CRS for union
  if (sf::st_crs(source_mask) != sf::st_crs(target_mask)) {
    target_mask <- sf::st_transform(target_mask, sf::st_crs(source_mask))
  }

  # Repair invalid geometries (common after raster->polygon + simplify)
  source_mask <- sf::st_make_valid(source_mask)
  target_mask <- sf::st_make_valid(target_mask)

  # Drop non-polygonal leftovers if any appear after repair
  source_mask <- cfcore_as_polygon(source_mask)
  target_mask <- cfcore_as_polygon(target_mask)

  # Union can still fail if there are tiny slivers; a zero-width buffer often fixes it
  union_sf <- tryCatch(
    sf::st_union(source_mask, target_mask),
    error = function(e) {
      sm <- sf::st_buffer(source_mask, 0)
      tm <- sf::st_buffer(target_mask, 0)
      sf::st_union(sm, tm)
    }
  )

  # Convert to terra mask
  union_vect <- terra::vect(union_sf)

  # Crop to union extent
  source <- terra::crop(source, terra::ext(union_vect))
  target <- terra::crop(target, terra::ext(union_vect))

  # Resample source onto target grid
  source <- terra::resample(source, target, method = method)

  # Mask both
  source <- terra::mask(source, union_vect)
  target <- terra::mask(target, union_vect)

  list(source = source, target = target, vect_mask = union_vect)
}


# ---- Continuous differences --------------------------------------------------

#' Compute a continuous difference raster
#'
#' Computes `later - earlier` and returns a continuous `terra::SpatRaster`
#' of change values (same units as the input rasters).
#'
#' @param earlier `terra::SpatRaster` (earlier epoch).
#' @param later `terra::SpatRaster` (later epoch).
#' @return `terra::SpatRaster` of continuous difference values.
#' @keywords internal
diff_continuous <- function(earlier, later) {
  if (!inherits(earlier, "SpatRaster") || !inherits(later, "SpatRaster")) {
    stop("`earlier` and `later` must be terra::SpatRaster.", call. = FALSE)
  }
  later - earlier
}

#' Continuous nDSM difference (later - earlier)
#'
#' @param earlier `terra::SpatRaster` (earlier epoch).
#' @param later `terra::SpatRaster` (later epoch).
#' @return `terra::SpatRaster` of continuous change values.
#' @keywords internal
diff_ndsm <- function(earlier, later) {
  diff_continuous(earlier, later)
}

#' Continuous DTM difference (later - earlier)
#'
#' @param earlier `terra::SpatRaster` (earlier epoch).
#' @param later `terra::SpatRaster` (later epoch).
#' @return `terra::SpatRaster` of continuous change values.
#' @keywords internal
diff_dtm <- function(earlier, later) {
  diff_continuous(earlier, later)
}

# ---- Reporting summaries (optional) -----------------------------------------

cfcore_bin_labels <- function(breaks, unit = "m") {
  cfcore_validate_breaks(breaks)

  fmt <- function(x) base::format(x, trim = TRUE, scientific = FALSE)
  unit_str <- if (is.null(unit) || !nzchar(unit)) "" else base::paste0(" ", unit)

  labs <- character(base::length(breaks) + 1L)

  labs[1] <- base::paste0("< ", fmt(breaks[1]), unit_str)

  if (base::length(breaks) > 1L) {
    for (i in 2:base::length(breaks)) {
      labs[i] <- base::paste0(fmt(breaks[i - 1L]), " to ", fmt(breaks[i]), unit_str)
    }
  }

  labs[base::length(labs)] <- base::paste0("> ", fmt(breaks[base::length(breaks)]), unit_str)
  labs
}


#' Validate breaks for binning
#'
#' @param breaks Numeric vector of strictly increasing finite breakpoints.
#' @keywords internal
cfcore_validate_breaks <- function(breaks) {
  if (!is.numeric(breaks) || base::length(breaks) < 1L) {
    stop("`breaks` must be a numeric vector with at least one element.", call. = FALSE)
  }
  if (base::any(!is.finite(breaks))) {
    stop("`breaks` must contain only finite values.", call. = FALSE)
  }
  if (base::is.unsorted(breaks, strictly = TRUE)) {
    stop("`breaks` must be strictly increasing.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Summarize a continuous difference raster by breaks
#'
#' Produces counts and percentages per bin defined by `breaks` (reporting).
#'
#' @param diff_raster `terra::SpatRaster` of continuous change values.
#' @param breaks Numeric vector of strictly increasing breakpoints.
#' @return A data.frame with `bin`, `count`, `percentage`.
#' @keywords internal
diff_summary_by_breaks <- function(diff_raster, breaks, unit = "m") {
  if (!inherits(diff_raster, "SpatRaster")) {
    stop("`diff_raster` must be terra::SpatRaster.", call. = FALSE)
  }
  cfcore_validate_breaks(breaks)

  v <- terra::values(diff_raster)
  v <- v[!base::is.na(v)]
  if (base::length(v) == 0L) {
    return(base::data.frame(bin = character(0), count = integer(0), percentage = numeric(0)))
  }

  edges <- c(-Inf, breaks, Inf)
  bins <- base::cut(v, breaks = edges, include.lowest = TRUE, right = TRUE)

  counts <- base::table(bins)
  pct <- (counts / base::sum(counts)) * 100

  labels <- cfcore_bin_labels(breaks, unit = unit)

  # Ensure ordering matches the bins
  out <- base::data.frame(
    bin = labels,
    count = base::as.integer(counts),
    percentage = base::as.numeric(pct),
    row.names = NULL
  )

  out
}

#' Continuous summary statistics for a difference raster
#'
#' Useful for QA/QC and API responses without binning.
#'
#' @param diff_raster `terra::SpatRaster` of continuous change values.
#' @param probs Numeric vector of quantile probabilities (default includes tails).
#' @return A list with `n`, `min`, `max`, `mean`, `sd`, `quantiles`.
#' @keywords internal
diff_summary_stats <- function(diff_raster, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)) {
  if (!inherits(diff_raster, "SpatRaster")) {
    stop("`diff_raster` must be terra::SpatRaster.", call. = FALSE)
  }

  v <- terra::values(diff_raster)
  v <- v[!base::is.na(v)]
  if (base::length(v) == 0L) {
    return(list(n = 0L, min = NA_real_, max = NA_real_, mean = NA_real_, sd = NA_real_, quantiles = numeric(0)))
  }

  list(
    n = base::length(v),
    min = base::min(v),
    max = base::max(v),
    mean = base::mean(v),
    sd = stats::sd(v),
    quantiles = stats::quantile(v, probs = probs, names = TRUE, na.rm = TRUE)
  )
}

# ---- Plotting helpers (optional) --------------------------------------------

#' Plot binned statistics for a continuous difference raster
#'
#' @param diff_raster `terra::SpatRaster` continuous change.
#' @param breaks Numeric vector of breakpoints.
#' @param palette Optional vector of colors (length = number of bins).
#' @param title Plot title.
#' @param xlab X-axis label.
#' @return A `ggplot2` object.
#' @keywords internal
plot_diff_binned <- function(diff_raster, breaks,
                             palette = NULL,
                             title = "Difference Raster Statistics",
                             xlab = "Change bin") {
  df <- diff_summary_by_breaks(diff_raster, breaks)
  df$bin <- factor(df$bin, levels = df$bin)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = bin,
      y = percentage,
      fill = bin
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_brewer(
      palette = "Set2",
      name = "Change class"   # <- LEGEND TITLE FIX
    ) +
    ggplot2::labs(
      title = title,
      x = "Vertical change range (m)",   # <- X AXIS
      y = "Area proportion (%)"           # <- Y AXIS
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 0)
    )
}

#' Plot a histogram for a continuous difference raster
#'
#' This is the pure continuous view (no bins by semantic breaks).
#'
#' @param diff_raster `terra::SpatRaster` continuous change.
#' @param bins Integer number of histogram bins.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @return A `ggplot2` object.
#' @keywords internal
plot_diff_hist <- function(diff_raster, bins = 60,
                           title = "Difference Histogram",
                           xlab = "Change") {
  if (!inherits(diff_raster, "SpatRaster")) {
    stop("`diff_raster` must be terra::SpatRaster.", call. = FALSE)
  }
  v <- terra::values(diff_raster)
  v <- v[!base::is.na(v)]

  plot_data <- base::data.frame(value = v)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = bins, color = "black") +
    ggplot2::labs(x = xlab, y = "Count") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal(base_size = 15)
}

#' Check whether a SpatRaster is empty (all NA)
#'
#' @param raster `terra::SpatRaster`.
#' @return `TRUE` if all raster values are NA, otherwise `FALSE`.
#' @keywords internal
is_empty <- function(raster) {
  base::all(base::is.na(terra::values(raster)))
}
