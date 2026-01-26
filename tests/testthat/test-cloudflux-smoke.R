test_that("cloudFlux workflow runs on packaged TTP data", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  data("TTP_15", package = "CFCore")
  data("TTP_23", package = "CFCore")

  cc <- cloudFlux_new(
    source_las = TTP_15,
    target_las = TTP_23,
    epsg = 26917,
    resolution = 1
  )

  cc$run()

  expect_true(inherits(cc$ndsm_diff, "SpatRaster"))
  expect_true(inherits(cc$dtm_diff, "SpatRaster"))

  stats <- cc$summary_stats()
  expect_true(is.list(stats))
  expect_true(all(c("ndsm", "dtm") %in% names(stats)))
  expect_true(is.list(stats$ndsm))
  expect_true(is.list(stats$dtm))
})
