#' Downsampled LiDAR point cloud – Tommy Thompson Park (2015)
#'
#' A downsampled airborne LiDAR point cloud representing Tommy Thompson Park,
#' collected in 2015. The dataset is intended for examples, testing, and
#' demonstration of CloudFluxCore point-cloud workflows (masking, DTM/nDSM
#' derivation, and change detection).
#'
#' The data have been spatially clipped and aggressively downsampled to keep
#' package size small while preserving representative structure.
#'
#' @format A [`lidR::LAS`] object.
#' @source Internal example data derived from City of Toronto LiDAR.
#' @seealso [TTP_23]
#' @name TTP_15
#' @aliases TTP_15
#' @keywords datasets
NULL

"TTP_15"


#' Downsampled LiDAR point cloud – Tommy Thompson Park (2023)
#'
#' A downsampled airborne LiDAR point cloud representing Tommy Thompson Park,
#' collected in 2023. This dataset pairs with [TTP_15] to demonstrate
#' multi-temporal LiDAR change workflows in CloudFluxCore.
#'
#' @format A [`lidR::LAS`] object.
#' @source Internal example data derived from City of Toronto LiDAR.
#' @seealso [TTP_15]
#' @name TTP_23
#' @aliases TTP_23
#' @keywords datasets
NULL

"TTP_23"
