test_that("aggregate_depth_window_by_replicate() computes a thickness-weighted mean per replicate", {
  sim_long <- data.frame(
    mukey = "1", cokey = "1", simulation_number = 1,
    hzdept_r = c(0, 20), hzdepb_r = c(20, 50),
    db = c(1.2, 1.5)
  )
  # Window [0, 30]: horizon 1 fully inside (20 cm overlap), horizon 2 partially (10 cm overlap)
  result <- aggregate_depth_window_by_replicate(sim_long, 0, 30, "db")
  expected <- (1.2 * 20 + 1.5 * 10) / 30
  expect_equal(result$db, expected, tolerance = 1e-9)
  expect_equal(result$top, 0)
  expect_equal(result$bottom, 30)
})

test_that("aggregate_depth_window_by_replicate() excludes horizons with zero overlap", {
  sim_long <- data.frame(
    mukey = "1", cokey = "1", simulation_number = 1,
    hzdept_r = c(0, 50), hzdepb_r = c(20, 100),
    db = c(1.2, 1.8)
  )
  result <- aggregate_depth_window_by_replicate(sim_long, 0, 20, "db")
  expect_equal(result$db, 1.2, tolerance = 1e-9)
})

test_that("aggregate_depth_window_by_replicate() aggregates each replicate independently", {
  sim_long <- data.frame(
    mukey = "1", cokey = c("1", "2"), simulation_number = c(1, 1),
    hzdept_r = c(0, 0), hzdepb_r = c(20, 20),
    db = c(1.2, 1.6)
  )
  result <- aggregate_depth_window_by_replicate(sim_long, 0, 20, "db")
  expect_equal(nrow(result), 2)
  expect_setequal(result$db, c(1.2, 1.6))
})

test_that("property_to_sim_column() maps recognized property ids and errors on unknown ones", {
  expect_equal(property_to_sim_column("ph"), "ph")
  expect_equal(property_to_sim_column("bulk_density"), "db")
  expect_equal(property_to_sim_column("clay"), "clay_total")
  expect_error(property_to_sim_column("awc"), "no simulated column mapping")
})

test_that("maybe_adjust_soil_data_depth_trend() passes through cokeys with fewer than min_depths distinct depths", {
  sim_long <- data.frame(
    cokey = "1", hzdept_r = 0, hzdepb_r = 20, simulation_number = 1:5,
    db = stats::rnorm(5, 1.4, 0.05)
  )
  result <- maybe_adjust_soil_data_depth_trend(sim_long, "db", min_depths = 2)
  expect_equal(nrow(result), 5)
  expect_equal(sort(result$db), sort(sim_long$db))
})

test_that("maybe_adjust_soil_data_depth_trend() parallel = TRUE matches parallel = FALSE", {
  testthat::skip_if_not_installed("GPfit")
  set.seed(99)
  build_cokey <- function(id, depths) {
    data.frame(
      cokey = id, hzdept_r = depths, hzdepb_r = depths + 20,
      simulation_number = 1,
      db = 1.2 + 0.02 * depths + stats::rnorm(length(depths), sd = 0.02)
    )
  }
  sim_long <- rbind(
    build_cokey("1", c(0, 20, 50, 100)),
    build_cokey("2", c(0, 30, 60)),
    build_cokey("3", c(0, 20, 40, 80, 120))
  )

  res_seq <- maybe_adjust_soil_data_depth_trend(sim_long, "db", min_depths = 2, parallel = FALSE)
  res_par <- maybe_adjust_soil_data_depth_trend(sim_long, "db", min_depths = 2, parallel = TRUE, n_cores = 2)

  order_cols <- c("cokey", "hzdept_r")
  res_seq <- res_seq[do.call(order, res_seq[order_cols]), ]
  res_par <- res_par[do.call(order, res_par[order_cols]), ]

  expect_equal(nrow(res_par), nrow(res_seq))
  expect_equal(res_par$cokey, res_seq$cokey)
  expect_equal(res_par$hzdept_r, res_seq$hzdept_r)
  expect_equal(res_par$db, res_seq$db, tolerance = 1e-6)
})

test_that("maybe_adjust_soil_data_depth_trend()/adjust_one_cokey_depth_trend() thread config through to reach joint_copula (Phase 10)", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 10: before this fix, config had no path from
  # simulate_ssurgo_mapunit_draws()'s own top-level API down to apply_gp_depth_trends()'s
  # dispatch, so "joint_copula" was unreachable except by calling apply_local_gp_adjustments()
  # directly (Phase 6's fix only closed that one hop). This confirms the full chain
  # (maybe_adjust_soil_data_depth_trend() -> adjust_one_cokey_depth_trend() ->
  # apply_local_gp_adjustments()) now carries config end to end.
  testthat::skip_if_not_installed("GPfit")
  set.seed(101)
  # Two properties, not one - apply_gp_depth_trends()'s vertical-correlation method dispatch only
  # branches with >= 2 properties; with just one, both methods trivially fall through the same
  # "individual adjustment" path and would never actually diverge (a real gap this test's first
  # draft had, caught the same way an identical gap was caught in Phase 11's own test).
  # Multiple realizations per depth (not a single simulation_number = 1) - with only one
  # realization, quantile-based remapping is a no-op for EITHER algorithm (nothing to remap
  # against), so both methods degenerate to the same unchanged output regardless of which is
  # used - a real degenerate case this test's first draft hit, not a source bug (confirmed:
  # preserve_correlation_structure_joint() already degrades gracefully - warns and falls back to
  # an identity property-correlation matrix - exactly as designed for too-little-data cases).
  depths <- c(0, 20, 50, 100)
  n_sims <- 30
  sim_long <- data.frame(
    cokey = "1",
    hzdept_r = rep(depths, each = n_sims), hzdepb_r = rep(depths + 20, each = n_sims),
    simulation_number = rep(seq_len(n_sims), times = length(depths)),
    db = 1.2 + 0.02 * rep(depths, each = n_sims) + stats::rnorm(n_sims * length(depths), sd = 0.02),
    wr_3b = 0.25 - 0.001 * rep(depths, each = n_sims) + stats::rnorm(n_sims * length(depths), sd = 0.01)
  )
  properties <- c("db", "wr_3b")

  retrofit_config <- get_default_configuration("validation")
  retrofit_config$monte_carlo$vertical_correlation_method <- "gp_quantile_retrofit"
  joint_config <- get_default_configuration("validation")
  joint_config$monte_carlo$vertical_correlation_method <- "joint_copula"

  set.seed(102)
  result_retrofit <- maybe_adjust_soil_data_depth_trend(sim_long, properties, min_depths = 2, config = retrofit_config)
  set.seed(102)
  result_joint <- maybe_adjust_soil_data_depth_trend(sim_long, properties, min_depths = 2, config = joint_config)

  expect_equal(nrow(result_joint), nrow(sim_long))
  expect_false(isTRUE(all.equal(result_retrofit$db, result_joint$db)))

  # config = NULL (the default) now reproduces the joint_copula result - Phase 13 flipped the
  # default; "gp_quantile_retrofit" remains reachable only as an explicit opt-out (already
  # confirmed via `retrofit_config` above).
  set.seed(102)
  result_default <- maybe_adjust_soil_data_depth_trend(sim_long, properties, min_depths = 2)
  expect_identical(result_default, result_joint)
})

test_that("maybe_adjust_soil_data_depth_trend() warns and passes through when GPfit isn't installed", {
  testthat::skip_if(nzchar(system.file(package = "GPfit")), "GPfit is installed - guard path not exercised")
  sim_long <- data.frame(cokey = "1", hzdept_r = c(0, 20), hzdepb_r = c(20, 50), simulation_number = 1, db = c(1.2, 1.5))
  expect_warning(result <- maybe_adjust_soil_data_depth_trend(sim_long, "db"), "GPfit not installed")
  expect_equal(result, sim_long)
})

test_that("rasterize_mukey_percentiles() correctly attaches per-mukey percentile values to each cell", {
  mukey_raster <- terra::rast(nrows = 3, ncols = 3, vals = c(101, 101, 102, 101, 102, 102, 103, 103, 103))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  percentile_by_mukey <- data.frame(mukey = c(101, 102, 103), P50 = c(1.1, 2.2, 3.3))
  result <- rasterize_mukey_percentiles(mukey_raster, percentile_by_mukey)

  expect_named(result, "P50")
  vals <- terra::values(result$P50)[, 1]
  expect_equal(vals, c(1.1, 1.1, 2.2, 1.1, 2.2, 2.2, 3.3, 3.3, 3.3))
})

test_that("rasterize_mukey_percentiles() handles multiple percentile columns", {
  mukey_raster <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 1, 2, 2))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  percentile_by_mukey <- data.frame(mukey = c(1, 2), P05 = c(0.5, 1.5), P95 = c(2.5, 3.5))
  result <- rasterize_mukey_percentiles(mukey_raster, percentile_by_mukey)

  expect_setequal(names(result), c("P05", "P95"))
  expect_equal(terra::values(result$P05)[, 1], c(0.5, 0.5, 1.5, 1.5))
  expect_equal(terra::values(result$P95)[, 1], c(2.5, 2.5, 3.5, 3.5))
})

test_that("fetch_ssurgo_mukey_raster()/simulate_ssurgo_mapunit_draws()/fetch_ssurgo_percentiles() require the live SDA service", {
  testthat::skip_if_offline()
  testthat::skip("Live NRCS Soil Data Access queries are not exercised in automated tests - see test-ssurgo-acquisition.R for the established precedent.")
})
