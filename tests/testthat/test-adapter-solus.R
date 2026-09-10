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

# ---------------------------------------------------------------------------
# S1 - fetch_solus_low_pred_high_multi() / fetch_solus_percentiles_multi(): batched multi-
# variable fetchSOLUS() (offline, mocked soilDB::fetchSOLUS()).
# ---------------------------------------------------------------------------

# Builds a synthetic multi-layer SpatRaster following the confirmed live layer-naming convention
# (paste0(var, "_", depth, "_cm_", suffix)), one constant-value layer per (variable, depth,
# suffix) combination - mirrors the live-verified 2026-09-04 spike shape (nlyr = n_vars * n_depths
# * n_suffixes).
.solus_multi_stack <- function(vars, depths, suffixes = c("p", "l", "h"), drop = character(0)) {
  combos <- expand.grid(var = vars, depth = depths, suffix = suffixes, stringsAsFactors = FALSE)
  combos$name <- paste0(combos$var, "_", combos$depth, "_cm_", combos$suffix)
  combos <- combos[!combos$name %in% drop, ]
  layers <- lapply(seq_len(nrow(combos)), function(i) {
    v <- match(combos$suffix[i], c("l", "p", "h"))  # low < pred < high, deterministic
    terra::rast(nrows = 1, ncols = 1, vals = as.numeric(combos$depth[i]) + v)
  })
  r <- do.call(c, layers)
  names(r) <- combos$name
  r
}

test_that("fetch_solus_low_pred_high_multi() fetches every variable in one fetchSOLUS() call", {
  fetch_n <- 0L
  vars <- c("claytotal", "dbovendry")
  depths <- c(0, 5, 15)

  testthat::local_mocked_bindings(
    fetchSOLUS = function(aoi_vect, depth_slices, variables, output_type, grid = TRUE) {
      fetch_n <<- fetch_n + 1L
      expect_setequal(variables, vars)
      .solus_multi_stack(vars, depths)
    },
    .package = "soilDB"
  )

  res <- fetch_solus_low_pred_high_multi(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"),
                                         vars, top_depth = 0, bottom_depth = 15)

  expect_equal(fetch_n, 1L)
  expect_named(res, vars)
  expect_s4_class(res$claytotal$pred, "SpatRaster")
  expect_s4_class(res$dbovendry$low, "SpatRaster")
  expect_true(all(terra::values(res$claytotal$low) <= terra::values(res$claytotal$pred)))
  expect_true(all(terra::values(res$claytotal$pred) <= terra::values(res$claytotal$high)))
})

test_that("fetch_solus_low_pred_high_multi() flags only the affected (variable, output_type) on partial failure", {
  vars <- c("claytotal", "dbovendry")
  depths <- c(0, 5, 15)
  # Drop every "high" layer for dbovendry only - a partial failure confined to one variable/output_type.
  missing <- paste0("dbovendry_", depths, "_cm_h")

  testthat::local_mocked_bindings(
    fetchSOLUS = function(...) .solus_multi_stack(vars, depths, drop = missing),
    .package = "soilDB"
  )

  res <- suppressWarnings(fetch_solus_low_pred_high_multi(
    terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"), vars, top_depth = 0, bottom_depth = 15
  ))

  expect_s4_class(res$claytotal$pred, "SpatRaster")
  expect_s4_class(res$claytotal$high, "SpatRaster")
  expect_s4_class(res$dbovendry$pred, "SpatRaster")
  expect_null(res$dbovendry$high)  # only the affected combo is NULL
  expect_warning(
    fetch_solus_low_pred_high_multi(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"),
                                    vars, top_depth = 0, bottom_depth = 15),
    "dbovendry.*95% high"
  )
})

test_that("fetch_solus_low_pred_high_multi() returns all-NULL per variable (not a crash) when fetchSOLUS() itself errors", {
  vars <- c("claytotal", "dbovendry")

  testthat::local_mocked_bindings(
    fetchSOLUS = function(...) stop("simulated network failure"),
    .package = "soilDB"
  )

  res <- suppressWarnings(fetch_solus_low_pred_high_multi(
    terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"), vars, top_depth = 0, bottom_depth = 15
  ))

  expect_named(res, vars)
  expect_null(res$claytotal$pred)
  expect_null(res$dbovendry$pred)
  expect_warning(
    fetch_solus_low_pred_high_multi(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"),
                                    vars, top_depth = 0, bottom_depth = 15),
    "fetchSOLUS\\(\\) failed"
  )
})

test_that("fetch_solus_percentiles_multi() packages every variable's low/pred/high into values/probs, NULL per variable if incomplete", {
  vars <- c("claytotal", "dbovendry")
  depths <- c(0, 5, 15)
  missing <- paste0("dbovendry_", depths, "_cm_h")

  testthat::local_mocked_bindings(
    fetchSOLUS = function(...) .solus_multi_stack(vars, depths, drop = missing),
    .package = "soilDB"
  )

  res <- suppressWarnings(fetch_solus_percentiles_multi(
    terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"), vars, top_depth = 0, bottom_depth = 15
  ))

  expect_named(res, vars)
  expect_equal(res$claytotal$probs, c(0.025, 0.5, 0.975))
  expect_named(res$claytotal$values, c("P025", "P50", "P975"))
  expect_null(res$dbovendry)  # incomplete (missing "high") -> whole-variable NULL, matches fetch_solus_percentiles()
})

# ---------------------------------------------------------------------------
# S2 - fetch_solus_site_level() / fetch_solus_restriction_depth() (bedrock/restriction-depth
# truncation, MULTI_PROPERTY_FUSION_PLAN.md task S2). Offline, mocked soilDB::fetchSOLUS(),
# matching the S1 spike's confirmed site-level layer-naming convention
# ("<variable>_all_cm_<suffix>", one layer per output type, no depth-slice looping).
# ---------------------------------------------------------------------------

.solus_site_level_stack <- function(variable, values = c(p = 40, l = 30, h = 50),
                                    suffixes = c("p", "l", "h")) {
  layers <- lapply(suffixes, function(s) terra::rast(nrows = 1, ncols = 1, vals = values[[s]]))
  r <- do.call(c, layers)
  names(r) <- paste0(variable, "_all_cm_", suffixes)
  r
}

test_that("fetch_solus_site_level() issues one fetchSOLUS() call with depth_slices = 'all' and parses the site-level layer names", {
  fetch_n <- 0L
  testthat::local_mocked_bindings(
    fetchSOLUS = function(aoi_vect, depth_slices, variables, output_type, grid = TRUE) {
      fetch_n <<- fetch_n + 1L
      expect_equal(depth_slices, "all")
      expect_equal(variables, "anylithicdpt")
      .solus_site_level_stack("anylithicdpt")
    },
    .package = "soilDB"
  )

  res <- fetch_solus_site_level(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"), "anylithicdpt")

  expect_equal(fetch_n, 1L)
  expect_equal(terra::values(res$pred)[1], 40)
  expect_equal(terra::values(res$low)[1], 30)
  expect_equal(terra::values(res$high)[1], 50)
})

test_that("fetch_solus_site_level() respects a restricted output_types request", {
  testthat::local_mocked_bindings(
    fetchSOLUS = function(aoi_vect, depth_slices, variables, output_type, grid = TRUE) {
      expect_equal(output_type, "prediction")
      .solus_site_level_stack("anylithicdpt", suffixes = "p")
    },
    .package = "soilDB"
  )

  res <- fetch_solus_site_level(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"),
                                "anylithicdpt", output_types = "prediction")
  expect_equal(terra::values(res$pred)[1], 40)
  expect_null(res$low)
  expect_null(res$high)
})

test_that("fetch_solus_site_level() returns all-NULL (not a crash) when fetchSOLUS() itself errors", {
  testthat::local_mocked_bindings(
    fetchSOLUS = function(...) stop("simulated network failure"),
    .package = "soilDB"
  )
  res <- suppressWarnings(fetch_solus_site_level(
    terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"), "anylithicdpt"
  ))
  expect_null(res$pred); expect_null(res$low); expect_null(res$high)
})

test_that("fetch_solus_restriction_depth() defaults to anylithicdpt and accepts resdept as an alternative", {
  seen_variable <- NULL
  testthat::local_mocked_bindings(
    fetchSOLUS = function(aoi_vect, depth_slices, variables, output_type, grid = TRUE) {
      seen_variable <<- variables
      .solus_site_level_stack(variables, values = c(p = 40, l = 30, h = 50), suffixes = "p")
    },
    .package = "soilDB"
  )

  aoi <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  fetch_solus_restriction_depth(aoi)
  expect_equal(seen_variable, "anylithicdpt")

  fetch_solus_restriction_depth(aoi, variable = "resdept")
  expect_equal(seen_variable, "resdept")
})

test_that("fetch_solus_restriction_depth() recodes right-censored (>= 201cm) values to Inf, not a literal cutoff", {
  testthat::local_mocked_bindings(
    fetchSOLUS = function(aoi_vect, depth_slices, variables, output_type, grid = TRUE) {
      terra::rast(nrows = 1, ncols = 4, vals = c(50, 150, 201, 250),
                  names = paste0(variables, "_all_cm_p"))
    },
    .package = "soilDB"
  )

  res <- fetch_solus_restriction_depth(terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070"))
  vals <- terra::values(res)[, 1]
  expect_equal(vals[1:2], c(50, 150))          # below the censor point - passed through unchanged
  expect_true(all(is.infinite(vals[3:4])))     # >= SOLUS_RESTRICTION_CENSOR_CM - recoded to Inf
})

test_that("fetch_solus_restriction_depth() returns NULL (not a crash) when the underlying fetch fails", {
  testthat::local_mocked_bindings(
    fetchSOLUS = function(...) stop("simulated network failure"),
    .package = "soilDB"
  )
  res <- suppressWarnings(fetch_solus_restriction_depth(
    terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  ))
  expect_null(res)
})

test_that("SOLUS_RESTRICTION_CENSOR_CM is 201, per SOLUS100 documentation section 2.1.3", {
  expect_equal(SOLUS_RESTRICTION_CENSOR_CM, 201)
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
