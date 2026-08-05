test_that("closest_solus_depth_slice() snaps to the native depth nearest the window's midpoint", {
  expect_equal(closest_solus_depth_slice(0, 10), 5)   # midpoint 5 -> exactly on a native slice
  expect_equal(closest_solus_depth_slice(0, 20), 5)   # midpoint 10 -> closer to 5 than 15
  expect_equal(closest_solus_depth_slice(20, 40), 30) # midpoint 30 -> exact
  expect_equal(closest_solus_depth_slice(90, 110), 100)
})

test_that("closest_solus_depth_slice() respects a custom available_slices vector", {
  expect_equal(closest_solus_depth_slice(0, 10, available_slices = c(0, 100)), 0)
})

test_that("fetch_solus_low_pred_high()/fetch_solus_percentiles() require the live SOLUS service", {
  testthat::skip_if_offline()
  # Live-verified scenario (Salinas Valley, CA - confirmed real SOLUS100 coverage): claytotal at
  # the native 0 cm depth slice returns exactly the layer names this function assumes
  # (claytotal_0_cm_p/_l/_h). Sys.unsetenv("PROJ_LIB") is also done in .onLoad(), but a stray
  # PROJ_LIB pointing at an incompatible external GDAL/PROJ install (documented in
  # code_ref/reanalysis-platform/HANDOFF_NOTES.md) can poison terra's PROJ context before
  # soilSIM's own .onLoad() runs, if terra was already attached earlier in the same session (as
  # it always is under devtools::test()/R CMD check, which load every test file's package
  # dependencies up front) - redone here defensively right before constructing a real AOI.
  Sys.unsetenv("PROJ_LIB")
  aoi <- terra::vect(
    "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
    "AOI has no valid CRS in this session - likely the documented PROJ_LIB environment issue (see HANDOFF_NOTES.md); not a soilSIM defect."
  )

  # fetchSOLUS() hits a real, sometimes transiently-unavailable external grid-index service
  # (independently confirmed working via an isolated live spike outside the automated suite - see
  # [[project_soilsim_package]] memory) - tolerate a graceful NULL (this function's own documented
  # fallback contract on fetch failure) rather than asserting the live service always succeeds.
  lph <- suppressWarnings(fetch_solus_low_pred_high(aoi, "claytotal", top_depth = 0, bottom_depth = 5))
  testthat::skip_if(is.null(lph$pred), "fetchSOLUS() was unavailable in this run (live external service) - not a soilSIM defect.")

  expect_s4_class(lph$pred, "SpatRaster")
  expect_s4_class(lph$low, "SpatRaster")
  expect_s4_class(lph$high, "SpatRaster")
  expect_true(all(terra::values(lph$low) <= terra::values(lph$pred), na.rm = TRUE))
  expect_true(all(terra::values(lph$pred) <= terra::values(lph$high), na.rm = TRUE))

  result <- fetch_solus_percentiles(aoi, "claytotal", top_depth = 0, bottom_depth = 5)
  expect_equal(result$probs, c(0.025, 0.5, 0.975))
  expect_named(result$values, c("P025", "P50", "P975"))
})
