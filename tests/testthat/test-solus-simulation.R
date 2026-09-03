test_that("closest_solus_depth_slice() snaps to the native depth nearest the window's midpoint", {
  expect_equal(closest_solus_depth_slice(0, 10), 5)   # midpoint 5 -> exactly on a native slice
  expect_equal(closest_solus_depth_slice(0, 20), 5)   # midpoint 10 -> closer to 5 than 15
  expect_equal(closest_solus_depth_slice(20, 40), 30) # midpoint 30 -> exact
  expect_equal(closest_solus_depth_slice(90, 110), 100)
})

test_that("closest_solus_depth_slice() respects a custom available_slices vector", {
  expect_equal(closest_solus_depth_slice(0, 10, available_slices = c(0, 100)), 0)
})

# ---------------------------------------------------------------------------
# solus_depth_window_weights() - trapezoidal depth-average over native slices
# ---------------------------------------------------------------------------

test_that("solus_depth_window_weights() sums to 1 and is sparse outside the window", {
  w <- solus_depth_window_weights(0, 30)
  expect_equal(sum(w), 1)
  expect_named(w, c("0", "5", "15", "30", "60", "100", "150"))
  expect_equal(unname(w[c("60", "100", "150")]), c(0, 0, 0))
  # 0-30 trapezoid over knots 0,5,15,30 -> raw 2.5/7.5/12.5/7.5, normalized by 30
  expect_equal(unname(w[c("0", "5", "15", "30")]), c(2.5, 7.5, 12.5, 7.5) / 30)
})

test_that("solus_depth_window_weights() is exact for a linear depth profile (== midpoint value)", {
  slices <- c(0, 5, 15, 30, 60, 100, 150)
  f <- function(d) 12 + 0.4 * d          # arbitrary linear profile
  for (win in list(c(0, 10), c(0, 30), c(10, 40), c(5, 100), c(0, 150))) {
    w <- solus_depth_window_weights(win[1], win[2])
    got <- sum(w * f(slices))
    expect_equal(got, f(mean(win)), tolerance = 1e-9,
                 info = sprintf("window %d-%d", win[1], win[2]))
  }
})

test_that("solus_depth_window_weights() interpolates the window end between bracketing slices", {
  w <- solus_depth_window_weights(0, 10)  # end at 10 -> halfway between slices 5 and 15
  expect_equal(sum(w), 1)
  expect_equal(unname(w[c("0", "5", "15")]), c(2.5, 6.25, 1.25) / 10)
})

test_that("solus_depth_window_weights() narrower-than-slice-spacing window averages the two bracketing slices", {
  w <- solus_depth_window_weights(0, 5)
  expect_equal(unname(w[c("0", "5")]), c(0.5, 0.5))
})

test_that("solus_depth_window_weights() clamps a window past the deepest slice", {
  w <- solus_depth_window_weights(140, 220)  # span clamps to 140-150; avg ~ value at 145
  expect_equal(sum(w), 1)
  expect_equal(unname(w[c("100", "150")]), c(0.1, 0.9))
})

test_that("solus_depth_window_weights() falls back to the nearest slice for a window entirely below the deepest", {
  w <- solus_depth_window_weights(160, 200)
  expect_equal(sum(w), 1)
  expect_equal(unname(w["150"]), 1)
})

test_that("solus_depth_window_weights() rejects a degenerate window", {
  expect_error(solus_depth_window_weights(10, 10), "bottom_depth > top_depth")
  expect_error(solus_depth_window_weights(30, 10), "bottom_depth > top_depth")
})

test_that("fetch_solus_low_pred_high()/fetch_solus_percentiles() require the live SOLUS service", {
  testthat::skip_if_offline()
  # Live-verified scenario (Salinas Valley, CA - confirmed real SOLUS100 coverage): claytotal over
  # 0-5 cm now pulls the native 0 cm and 5 cm slices (claytotal_0_cm_p/_l/_h + claytotal_5_cm_*)
  # and depth-averages them via solus_depth_window_weights(). Sys.unsetenv("PROJ_LIB") is also done in .onLoad(), but a stray
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
