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
