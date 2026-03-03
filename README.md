---
CFCore
---

<!-- badges: start -->
[![R-CMD-check](https://github.com/cscarpon/CFCore/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cscarpon/CFCore/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

**CFCore** is the core analytical engine behind *CloudFlux*.  
It provides reusable, API-friendly R functions and reference classes for **LiDAR-based change detection**, including point-cloud preprocessing, raster generation, raster alignment, and continuous change analysis.

CFCore is intentionally **UI-agnostic**. Functions and classes are designed to be called directly from R scripts, batch pipelines, reproducible research workflows, and downstream services.

---

## Architecture

```mermaid
flowchart TD
  A[Inputs: LAS/LAZ] --> B[Mask + Denoise]
  B --> C[DTM + nDSM]
  C --> D[Align rasters]
  D --> E[Diff rasters]
  E --> F[Stats / Plots / Map]
```

## Key features

- LiDAR point-cloud preprocessing
  - Noise filtering (SOR-based)
  - 2D footprint / mask generation
- Raster derivation
  - DTM and nDSM generation
  - CRS handling and reprojection
- Two-epoch change detection
  - Raster alignment to a shared union mask
  - Continuous difference rasters (nDSM, DTM)
- Summary and reporting
  - Continuous summary statistics
  - Optional binned interpretation summaries
  - Histograms and binned plots
- Optional ICP alignment
  - Open3D-based ICP via Python (`reticulate`)
- High-level orchestration
  - `cloudFlux` reference class for end-to-end workflows

---

## Installation

### Install CFCore from GitHub

``` r
# install.packages("pak")
pak::pak("cscarpon/CFCore")
```

## installing Python environment for ICP alignment 
```
conda env create -f inst/py/conda-env.yml
conda activate icp_conda
```

## Example code and High level workflow with `cloudflux`

```r
library(CFCore)

# Load example datasets (Tommy Thompson Park 2015 & 2023)
data("TTP_15", package = "CFCore")
data("TTP_23", package = "CFCore")

# Initialize the CloudFlux workflow
cc <- cloudFlux_new(
  source_las = TTP_15,
  target_las = TTP_23,
  epsg = 26917,
  resolution = 1,
  use_icp = TRUE,     # Automatically run Open3D ICP alignment during execution
  use_gpu = TRUE      # Toggle CUDA GPU acceleration (falls back to CPU if unavailable)
)

# Execute the full workflow (Alignment -> Masking -> Rasterization -> Differencing)
cc$run()

# Continuous summary statistics
cc$summary_stats()

# Optional binned interpretation (e.g., minimal change, moderate growth/loss, large changes)
breaks <- c(-10, -0.5, 0.5, 10)
cc$summary_by_breaks(breaks)

# Visualizations
cc$plot_hist(which = "ndsm", title = "Canopy Height Change Distribution")
cc$plot_binned(breaks, which = "ndsm", title = "Binned Canopy Change")

# Interactive Leaflet map (Great for RStudio viewer or HTML Markdown)
m <- cc$map(diff_breaks = breaks)
m 

# Export the generated DTMs, nDSMs, and difference rasters to a local directory
cc$save_rasters(out_dir = "cfcore_outputs", prefix = "ttp")
```

```r
## Workflow with unique paths:

cc <- cloudFlux_new(
  source_path = "path/to/TTP_15.laz",
  target_path = "path/to/TTP_23.laz",
  epsg = 26917
)

cc$run()
```

```r

## Using the icp-alignment workflow:

cc <- cloudFlux_new(
  source_las = TTP_15,
  target_las = TTP_23,
  epsg = 26917,
  icp_align = TRUE,
  icp_condaenv = "icp_conda",
  voxel_size = 0.05,
  icp_method = "point-to-plane",
  icp_use_gpu = TRUE,   # optional (default TRUE)
  icp_max_iter = 50     # optional (default 50)
)

cc$run()
```
