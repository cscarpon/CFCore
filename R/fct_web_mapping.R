# ---- internal helpers --------------------------------------------------------

cfcore_epsg_string <- function(epsg) {
  base::paste0("EPSG:", base::as.integer(epsg))
}

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

cfcore_labels_from_breaks <- function(breaks) {
  cfcore_validate_breaks(breaks)
  edges <- c(-Inf, breaks, Inf)

  fmt <- function(x) {
    if (is.infinite(x)) {
      if (x < 0) "-Inf" else "Inf"
    } else {
      base::format(x, trim = TRUE, scientific = FALSE)
    }
  }

  labs <- character(base::length(edges) - 1L)
  for (i in base::seq_along(labs)) {
    a <- edges[i]
    b <- edges[i + 1L]
    labs[i] <- if (is.infinite(a)) {
      base::paste0("< ", fmt(b))
    } else if (is.infinite(b)) {
      base::paste0("> ", fmt(a))
    } else {
      base::paste0(fmt(a), " to ", fmt(b))
    }
  }
  labs
}

cfcore_mask_to_vect <- function(mask, epsg = 4326) {
  if (inherits(mask, "SpatVector")) {
    terra::project(mask, cfcore_epsg_string(epsg))
  } else if (inherits(mask, "sf")) {
    terra::vect(sf::st_transform(mask, epsg))
  } else {
    stop("mask must be a terra::SpatVector or sf object.", call. = FALSE)
  }
}

cfcore_mask_to_sf <- function(mask, epsg = 4326) {
  if (inherits(mask, "SpatVector")) {
    sf::st_as_sf(terra::project(mask, cfcore_epsg_string(epsg)))
  } else if (inherits(mask, "sf")) {
    sf::st_transform(mask, epsg)
  } else {
    stop("mask must be a terra::SpatVector or sf object.", call. = FALSE)
  }
}

cfcore_project_mask <- function(mask, epsg = 4326) {
  cfcore_mask_to_vect(mask, epsg = epsg)
}

cfcore_values_non_na <- function(raster) {
  v <- terra::values(raster)
  v[!base::is.na(v)]
}

cfcore_diff_domain <- function(diff_raster, diff_range = NULL, diff_quantiles = c(0.02, 0.98)) {
  v <- cfcore_values_non_na(diff_raster)
  if (base::length(v) == 0L) return(NULL)

  if (!is.null(diff_range)) {
    if (!is.numeric(diff_range) || base::length(diff_range) != 2L) {
      stop("`diff_range` must be numeric length-2 c(min, max).", call. = FALSE)
    }
    return(diff_range)
  }

  if (!is.numeric(diff_quantiles) || base::length(diff_quantiles) != 2L) {
    stop("`diff_quantiles` must be numeric length-2, e.g. c(0.02, 0.98).", call. = FALSE)
  }

  stats::quantile(v, probs = diff_quantiles, na.rm = TRUE, names = FALSE)
}

# ---- mapping functions -------------------------------------------------------

#' Display a leaflet map with rasters and optional continuous change layers
#'
#' Builds a leaflet map containing DTM and nDSM rasters for two epochs, optional
#' continuous difference rasters (`later - earlier`), and an optional AOI mask.
#'
#' @param dtm1,ndsm1,dtm2,ndsm2 `terra::SpatRaster` layers (may be `NULL`).
#' @param dtm_diff,ndsm_diff Optional continuous difference rasters (may be `NULL`).
#' @param area_mask Optional AOI mask. Prefer `terra::SpatVector` or `sf`.
#' @param epsg EPSG integer for web display. Defaults to 4326.
#' @param diff_palette Palette name or vector accepted by `leaflet::colorNumeric()`.
#'   If `NULL`, uses `"viridis"`.
#' @param diff_range Optional numeric length-2 c(min, max) for diff legend/rendering.
#'   If `NULL`, `diff_quantiles` are used.
#' @param diff_quantiles Numeric length-2 quantiles used to set diff rendering domain
#'   when `diff_range` is `NULL`. Defaults to c(0.02, 0.98) to reduce outlier impact.
#' @param diff_breaks Optional numeric vector of breakpoints for an additional binned
#'   legend (labels only). Does not alter raster rendering.
#' @param diff_break_colors Optional colors for the binned legend (length = bins).
#'   If `NULL`, a default qualitative palette is generated.
#' @param measure_colors Length-2 character vector for measure tool colors:
#'   `c(activeColor, completedColor)`.
#' @return A `leaflet` map widget.
#' @keywords internal
displayMap <- function(dtm1, ndsm1, dtm2, ndsm2, dtm_diff, ndsm_diff, area_mask,
                       epsg = 4326,
                       diff_palette = NULL,
                       diff_range = NULL,
                       diff_quantiles = c(0.02, 0.98),
                       diff_breaks = NULL,
                       diff_break_colors = NULL,
                       measure_colors = c("#3D535D", "#7D4479")) {
  m <- leaflet::leaflet() |>
    leaflet::addTiles() |>
    leaflet::addScaleBar(position = "bottomleft") |>
    leaflet::addMeasure(
      position = "topleft",
      primaryLengthUnit = "meters",
      primaryAreaUnit = "sqmeters",
      activeColor = measure_colors[1],
      completedColor = measure_colors[2]
    )

  if (!is.null(area_mask)) {
    area_mask <- cfcore_project_mask(area_mask, epsg = epsg)
  }

  dtm1_m <- processRasterMask(dtm1, area_mask, epsg = epsg)
  ndsm1_m <- processRasterMask(ndsm1, area_mask, epsg = epsg)
  dtm2_m  <- processRasterMask(dtm2, area_mask, epsg = epsg)
  ndsm2_m <- processRasterMask(ndsm2, area_mask, epsg = epsg)

  ndsm_diff_m <- processRasterMask(ndsm_diff, area_mask, epsg = epsg)
  dtm_diff_m  <- processRasterMask(dtm_diff,  area_mask, epsg = epsg)

  # Rasters
  m <- addRasterLayer(m, dtm1_m,  "DTM (Source)",  "magma")
  m <- addRasterLayer(m, dtm2_m,  "DTM (Target)",  "magma")
  m <- addRasterLayer(m, ndsm1_m, "nDSM (Source)", "magma")
  m <- addRasterLayer(m, ndsm2_m, "nDSM (Target)", "magma")

  if (!is.null(ndsm_diff_m)) {
    m <- addDiffRasterLayer(
      m, ndsm_diff_m,
      group = "Difference nDSM",
      palette = diff_palette,
      diff_range = diff_range,
      diff_quantiles = diff_quantiles
    )
  }
  if (!is.null(dtm_diff_m)) {
    m <- addDiffRasterLayer(
      m, dtm_diff_m,
      group = "Difference DTM",
      palette = diff_palette,
      diff_range = diff_range,
      diff_quantiles = diff_quantiles
    )
  }

  if (!is.null(area_mask)) {
    m <- leaflet::addPolygons(
      m,
      data = cfcore_mask_to_sf(area_mask, epsg = epsg),
      color = "black",
      fill = FALSE,
      group = "Mask"
    )
  }

  overlayGroups <- character(0)
  if (!is.null(dtm1_m))      overlayGroups <- c(overlayGroups, "DTM (Source)")
  if (!is.null(dtm2_m))      overlayGroups <- c(overlayGroups, "DTM (Target)")
  if (!is.null(ndsm1_m))     overlayGroups <- c(overlayGroups, "nDSM (Source)")
  if (!is.null(ndsm2_m))     overlayGroups <- c(overlayGroups, "nDSM (Target)")
  if (!is.null(dtm_diff_m))  overlayGroups <- c(overlayGroups, "Difference DTM")
  if (!is.null(ndsm_diff_m)) overlayGroups <- c(overlayGroups, "Difference nDSM")
  if (!is.null(area_mask))   overlayGroups <- c(overlayGroups, "Mask")

  m <- leaflet::addLayersControl(
    m,
    overlayGroups = overlayGroups,
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )

  # Default visibility: only Mask on
  for (g in overlayGroups) {
    if (!identical(g, "Mask")) m <- leaflet::hideGroup(m, g)
  }
  if ("Mask" %in% overlayGroups) m <- leaflet::showGroup(m, "Mask")

  # Legends: create all, then JS will manage visibility.
  if (!is.null(dtm1_m))  m <- addLegendLayer(m, dtm1_m,  "DTM (Source)",  "dtm1Legend",  "magma")
  if (!is.null(dtm2_m))  m <- addLegendLayer(m, dtm2_m,  "DTM (Target)",  "dtm2Legend",  "magma")
  if (!is.null(ndsm1_m)) m <- addLegendLayer(m, ndsm1_m, "nDSM (Source)", "ndsm1Legend", "magma")
  if (!is.null(ndsm2_m)) m <- addLegendLayer(m, ndsm2_m, "nDSM (Target)", "ndsm2Legend", "magma")

  # IMPORTANT: do NOT pass `group=` to addDiffLegend() unless its signature supports it.
  if (!is.null(ndsm_diff_m)) {
    m <- addDiffLegend(
      m, ndsm_diff_m,
      title = "Change in Normalized Height\nModel (m)",
      layerId = "diff_ndsm_legend",
      palette = diff_palette,
      diff_range = diff_range,
      diff_quantiles = diff_quantiles
    )
  }
  if (!is.null(dtm_diff_m)) {
    m <- addDiffLegend(
      m, dtm_diff_m,
      title = "Change in Digital Terrain\nModel (m)",
      layerId = "diff_dtm_legend",
      palette = diff_palette,
      diff_range = diff_range,
      diff_quantiles = diff_quantiles
    )
  }

  if (!is.null(diff_breaks)) {
    cfcore_validate_breaks(diff_breaks)
    labels <- cfcore_labels_from_breaks(diff_breaks)

    k <- base::length(labels)
    colors <- diff_break_colors
    if (is.null(colors)) {
      colors <- grDevices::hcl.colors(k, palette = "Set 2")
    } else if (base::length(colors) != k) {
      stop("`diff_break_colors` length must match number of bins: ", k, call. = FALSE)
    }

    m <- leaflet::addLegend(
      m,
      colors = colors,
      labels = labels,
      position = "bottomleft",
      title = "Interpretation bins",
      layerId = "diffBinsLegend",
      opacity = 1
    )
  }

  # JS: assign stable DOM ids to legend controls by matching title text, then toggle.
  # NOTE: JS normalizes whitespace, so "\n" becomes a space for matching.
  legend_titles_to_ids <- list(
    "DTM (Source)" = "dtm1Legend",
    "DTM (Target)" = "dtm2Legend",
    "nDSM (Source)" = "ndsm1Legend",
    "nDSM (Target)" = "ndsm2Legend",
    "Change in Normalized Height Model (m)" = "diff_ndsm_legend",
    "Change in Digital Terrain Model (m)" = "diff_dtm_legend",
    "Interpretation bins" = "diffBinsLegend"
  )

  overlay_to_legend_id <- list(
    "DTM (Source)" = "dtm1Legend",
    "DTM (Target)" = "dtm2Legend",
    "nDSM (Source)" = "ndsm1Legend",
    "nDSM (Target)" = "ndsm2Legend",
    "Difference nDSM" = "diff_ndsm_legend",
    "Difference DTM" = "diff_dtm_legend"
  )

  m <- htmlwidgets::onRender(
    m,
    "
    function(el, x, data) {
      var map = this;

      function legendNodes() {
        return Array.prototype.slice.call(
          el.querySelectorAll('.leaflet-control.info.legend')
        );
      }

      function normalize(s) {
        return String(s || '').replace(/\\s+/g, ' ').trim();
      }

      function assignLegendIdsByTitle() {
        var titleToId = data.legend_titles_to_ids || {};
        legendNodes().forEach(function(node) {
          if (node.id) return;
          var txt = normalize(node.textContent);

          Object.keys(titleToId).forEach(function(title) {
            if (node.id) return;
            var t = normalize(title);
            if (txt.indexOf(t) !== -1) node.id = titleToId[title];
          });
        });
      }

      function hideAllLegends() {
        (data.all_legend_ids || []).forEach(function(id) {
          var n = document.getElementById(id);
          if (n) n.style.display = 'none';
        });
      }

      function showLegend(id) {
        if (!id) return;
        var n = document.getElementById(id);
        if (n) n.style.display = '';
      }

      function groupIsActive(groupName) {
        if (!map.layerManager || !map.layerManager._byGroup) return false;
        var grp = map.layerManager._byGroup[groupName];
        if (!grp) return false;
        var layerIds = Object.keys(grp);
        for (var i = 0; i < layerIds.length; i++) {
          var layer = grp[layerIds[i]];
          if (layer && map.hasLayer(layer)) return true;
        }
        return false;
      }

      function chooseActiveOverlay() {
        var overlayToLegend = data.overlay_to_legend_id || {};
        if (map._cf_lastOverlay && groupIsActive(map._cf_lastOverlay)) {
          return map._cf_lastOverlay;
        }
        var names = Object.keys(overlayToLegend);
        for (var i = 0; i < names.length; i++) {
          if (groupIsActive(names[i])) return names[i];
        }
        return null;
      }

      function updateLegends() {
        assignLegendIdsByTitle();
        hideAllLegends();

        var overlayToLegend = data.overlay_to_legend_id || {};
        var activeOverlay = chooseActiveOverlay();
        if (activeOverlay && overlayToLegend[activeOverlay]) {
          showLegend(overlayToLegend[activeOverlay]);
          if ((activeOverlay === 'Difference nDSM' || activeOverlay === 'Difference DTM') &&
              data.has_bins_legend) {
            showLegend('diffBinsLegend');
          }
        }
      }

      map.on('overlayadd', function(e) {
        map._cf_lastOverlay = e.name;
        updateLegends();
      });

      map.on('overlayremove', function(e) {
        updateLegends();
      });

      updateLegends();
      hideAllLegends();
    }
    ",
    data = list(
      legend_titles_to_ids = legend_titles_to_ids,
      overlay_to_legend_id = overlay_to_legend_id,
      all_legend_ids = unique(unname(unlist(legend_titles_to_ids))),
      has_bins_legend = !is.null(diff_breaks)
    )
  )

  m
}

#' Project and optionally mask a raster for web display
#'
#' @param raster `terra::SpatRaster` or `NULL`.
#' @param mask Optional AOI mask (`terra::SpatVector` or `sf`), projected to display CRS.
#' @param epsg EPSG integer for display CRS (defaults to 4326).
#' @return Projected (and masked) `terra::SpatRaster`, or `NULL`.
#' @keywords internal
processRasterMask <- function(raster, mask, epsg = 4326) {
  if (is.null(raster)) return(NULL)

  raster_m <- terra::project(raster, cfcore_epsg_string(epsg))

  if (!is.null(mask)) {
    raster_m <- terra::mask(raster_m, cfcore_mask_to_vect(mask, epsg = epsg))
  }

  raster_m
}

#' Add a continuous raster layer to a leaflet map
#'
#' @param map Leaflet map.
#' @param raster `terra::SpatRaster` or `NULL`.
#' @param name Group name for leaflet layer control.
#' @param color_palette Palette name or vector accepted by `leaflet::colorNumeric()`.
#' @return Leaflet map.
#' @keywords internal
addRasterLayer <- function(map, raster, name, color_palette) {
  if (is.null(raster)) return(map)

  vals <- cfcore_values_non_na(raster)
  if (base::length(vals) == 0L) return(map)

  pal <- leaflet::colorNumeric(color_palette, domain = vals, na.color = "transparent")

  leaflet::addRasterImage(map, raster, colors = pal, group = name, maxBytes = Inf, opacity = 1)
}

#' Add a legend for a continuous raster layer
#'
#' @param map Leaflet map.
#' @param raster `terra::SpatRaster` or `NULL`.
#' @param title Legend title.
#' @param layerId Leaflet layer ID.
#' @param color_palette Palette name or vector accepted by `leaflet::colorNumeric()`.
#' @return Leaflet map.
#' @keywords internal
addLegendLayer <- function(map, raster, title, layerId, color_palette) {
  if (is.null(raster)) return(map)

  vals <- cfcore_values_non_na(raster)
  if (base::length(vals) == 0L) return(map)

  pal <- leaflet::colorNumeric(color_palette, domain = vals, na.color = "transparent")

  leaflet::addLegend(
    map,
    pal = pal,
    values = vals,
    position = "bottomright",
    title = title,
    layerId = layerId,
    opacity = 1
  )
}

#' Add a continuous difference raster layer (with controlled domain)
#'
#' Uses `diff_range` if provided; otherwise uses quantiles to set the color scale
#' domain to reduce outlier impact.
#'
#' @param map Leaflet map.
#' @param diff_raster `terra::SpatRaster` continuous difference raster.
#' @param group Leaflet layer group name.
#' @param palette Palette name or vector accepted by `leaflet::colorNumeric()`.
#' @param diff_range Optional numeric length-2 c(min, max) for domain.
#' @param diff_quantiles Quantiles used for domain if `diff_range` is NULL.
#' @return Leaflet map.
#' @keywords internal
addDiffRasterLayer <- function(map, diff_raster, group,
                               palette = NULL,
                               diff_range = NULL,
                               diff_quantiles = c(0.02, 0.98)) {
  if (is.null(diff_raster)) return(map)

  if (is.null(palette)) palette <- "viridis"

  domain <- cfcore_diff_domain(diff_raster, diff_range = diff_range, diff_quantiles = diff_quantiles)
  if (is.null(domain)) return(map)

  pal <- leaflet::colorNumeric(palette = palette, domain = domain, na.color = "transparent")

  leaflet::addRasterImage(map, diff_raster, colors = pal, group = group, maxBytes = Inf, opacity = 1)
}

#' Add a continuous legend for a difference raster (with controlled domain)
#'
#' @param map Leaflet map.
#' @param diff_raster `terra::SpatRaster` continuous difference raster.
#' @param title Legend title.
#' @param layerId Leaflet layer ID.
#' @param palette Palette name or vector accepted by `leaflet::colorNumeric()`.
#' @param diff_range Optional numeric length-2 c(min, max) for domain.
#' @param diff_quantiles Quantiles used for domain if `diff_range` is NULL.
#' @return Leaflet map.
#' @keywords internal
addDiffLegend <- function(map, diff_raster, title, layerId,
                          palette = NULL,
                          diff_range = NULL,
                          diff_quantiles = c(0.02, 0.98)) {
  if (is.null(diff_raster)) return(map)
  if (is.null(palette)) palette <- "viridis"

  domain <- cfcore_diff_domain(diff_raster, diff_range = diff_range, diff_quantiles = diff_quantiles)
  if (is.null(domain)) return(map)

  pal <- leaflet::colorNumeric(palette = palette, domain = domain, na.color = "transparent")

  leaflet::addLegend(
    map,
    pal = pal,
    values = domain,
    position = "bottomright",
    title = title,
    layerId = layerId,
    opacity = 1
  )
}


#' Display index polygons for loaded spatial containers
#'
#' @param index `sf::sfc` or `sf` geometry representing an extent polygon.
#' @param epsg EPSG integer for display CRS (defaults to 4326).
#' @param measure_colors Length-2 character vector for measure tool colors.
#' @return A `leaflet` map widget.
#' @keywords internal
displayIndex <- function(index, epsg = 4326, measure_colors = c("#3D535D", "#7D4479")) {
  m <- leaflet::leaflet() |>
    leaflet::addTiles() |>
    leaflet::addScaleBar(position = "bottomleft") |>
    leaflet::addMeasure(
      position = "topleft",
      primaryLengthUnit = "meters",
      primaryAreaUnit = "sqmeters",
      activeColor = measure_colors[1],
      completedColor = measure_colors[2]
    )

  if (!is.null(index)) {
    index <- sf::st_transform(index, epsg)
    m <- leaflet::addPolygons(m, data = index, color = "black", fill = FALSE, group = "Index")
    m <- leaflet::addLayersControl(
      m,
      overlayGroups = "Index",
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )
    m <- leaflet::showGroup(m, "Index")
  }

  m
}
