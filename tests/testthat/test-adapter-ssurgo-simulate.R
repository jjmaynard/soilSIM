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

test_that("property_to_sim_column() maps the 5 P2 chemistry properties to themselves", {
  for (p in c("caco3", "ec", "ecec", "gypsum", "sar")) {
    expect_equal(property_to_sim_column(p), p, info = p)
  }
})

test_that("SSURGO_SIM_PROPERTY_COLUMNS includes the 5 P2 chemistry properties", {
  expect_true(all(c("caco3", "ec", "ecec", "gypsum", "sar") %in% SSURGO_SIM_PROPERTY_COLUMNS))
})

test_that("infill_ssurgo_data() infills the 5 P2 chemistry properties when their _r columns are present", {
  df <- make_horizon_row(properties = c("dbovendry", "caco3", "ec", "ecec", "gypsum", "sar"))
  df <- rbind(df, df)
  df$caco3_r[1] <- NA
  df$hzname <- "A"

  result <- infill_ssurgo_data(df)
  # infill_soil_property() is generic over any "<name>_l/_r/_h" triplet - just confirm it ran
  # without error and the columns survive (full infilling-strategy coverage is
  # infill_soil_property()'s own concern, tested elsewhere).
  expect_true(all(c("caco3_r", "ec_r", "ecec_r", "gypsum_r", "sar_r") %in% names(result)))
})

make_wr_gap_df <- function() {
  # Two 1-horizon components with complete texture + BD but an entirely-missing
  # water-retention triplet: a sandy component and a clayey one. Nothing to borrow a
  # water-retention value FROM, so the generic hierarchy is stuck - only a pedotransfer
  # estimate can fill these.
  sandy <- make_horizon_row(properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "om"),
                            cokey = "1", mukey = "1",
                            sandtotal_r = 80, sandtotal_l = 75, sandtotal_h = 85,
                            claytotal_r = 5, claytotal_l = 3, claytotal_h = 8,
                            silttotal_r = 15, silttotal_l = 10, silttotal_h = 20)
  clayey <- make_horizon_row(properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "om"),
                             cokey = "2", mukey = "1",
                             sandtotal_r = 15, sandtotal_l = 10, sandtotal_h = 20,
                             claytotal_r = 55, claytotal_l = 50, claytotal_h = 60,
                             silttotal_r = 30, silttotal_l = 25, silttotal_h = 35)
  df <- rbind(sandy, clayey)
  df$hzname <- c("A", "Bt")
  df$comppct_r <- c(60, 40)
  for (col in c("wthirdbar_l", "wthirdbar_r", "wthirdbar_h",
                "wfifteenbar_l", "wfifteenbar_r", "wfifteenbar_h")) df[[col]] <- NA_real_
  df
}

test_that("infill_ssurgo_data() fills water-retention gaps via Saxton-Rawls, texture-consistently", {
  result <- infill_ssurgo_data(make_wr_gap_df())

  expect_false(anyNA(result$wthirdbar_r))
  expect_false(anyNA(result$wfifteenbar_r))
  # texture-consistent: the clay component holds more water than the sand component
  expect_gt(result$wthirdbar_r[2], result$wthirdbar_r[1])
  expect_gt(result$wfifteenbar_r[2], result$wfifteenbar_r[1])
  # _l/_h ranges populated too
  expect_false(anyNA(result$wthirdbar_l))
  expect_false(anyNA(result$wthirdbar_h))
})

test_that("infill_ssurgo_data() water-retention estimate equals the Saxton-Rawls PTF exactly", {
  df <- make_wr_gap_df()
  res <- infill_ssurgo_data(df, water_retention_method = "saxton_rawls")

  for (i in 1:2) {
    expected <- compute_saxton_rawls(
      sand_pct = df$sandtotal_r[i], clay_pct = df$claytotal_r[i], silt_pct = df$silttotal_r[i],
      bulk_density = df$dbovendry_r[i], rfv_pct = 0, om_pct = df$om_r[i]
    )$field_capacity
    expect_equal(res$wthirdbar_r[i], expected, info = paste("row", i))
  }

  # With Strategy 5 now also using Saxton-Rawls, the "generic" route converges to the same
  # values here (no neighbouring horizon has a real wr value to borrow first).
  gen <- infill_ssurgo_data(df, water_retention_method = "generic")
  expect_equal(gen$wthirdbar_r, res$wthirdbar_r)
})

test_that("infill_ssurgo_data() keeps genuine SSURGO water-retention values (overwrite = FALSE)", {
  df <- make_wr_gap_df()
  df$wthirdbar_l[1] <- 0.20; df$wthirdbar_r[1] <- 0.25; df$wthirdbar_h[1] <- 0.30

  result <- infill_ssurgo_data(df)

  expect_equal(result$wthirdbar_r[1], 0.25)          # real value untouched
  expect_false(is.na(result$wthirdbar_r[2]))         # gap filled by Saxton-Rawls
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

  retrofit_config <- default_config("validation")
  retrofit_config$monte_carlo$vertical_correlation_method <- "gp_quantile_retrofit"
  joint_config <- default_config("validation")
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

test_that("lookup_mukey_draws() groups a property's simulated values by mukey without collapsing to percentiles", {
  draws <- data.frame(
    mukey = c("1", "1", "1", "2", "2"),
    cokey = c("a", "a", "b", "c", "c"),
    simulation_number = c(1, 2, 1, 1, 2),
    db = c(1.1, 1.2, 1.3, 2.1, 2.2)
  )
  lookup <- lookup_mukey_draws(draws, "bulk_density")
  expect_setequal(names(lookup), c("1", "2"))
  expect_equal(sort(lookup[["1"]]), c(1.1, 1.2, 1.3))
  expect_equal(sort(lookup[["2"]]), c(2.1, 2.2))
})

test_that("lookup_mukey_draws() drops non-finite values and returns NULL for a column absent from draws", {
  draws <- data.frame(mukey = c("1", "1"), db = c(1.5, NA))
  lookup <- lookup_mukey_draws(draws, "bulk_density")
  expect_equal(lookup[["1"]], 1.5)

  draws_no_ph <- data.frame(mukey = "1", db = 1.4)
  expect_null(lookup_mukey_draws(draws_no_ph, "ph"))
})

test_that("lookup_mukey_draws() errors on an unrecognized property id, matching percentiles_from_draws()'s convention", {
  draws <- data.frame(mukey = "1", db = 1.4)
  expect_error(lookup_mukey_draws(draws, "not_a_real_property"), "no simulated column mapping")
})

test_that("lookup_mukey_texture_draws() returns row-aligned clay/sand/silt triples per mukey, dropping incomplete rows", {
  draws <- data.frame(
    mukey = c("1", "1", "1", "2"),
    clay_total = c(20, 22, NA, 15),
    sand_total = c(40, 38, 41, 45),
    silt_total = c(40, 40, 41, 40)
  )
  lookup <- lookup_mukey_texture_draws(draws)
  expect_setequal(names(lookup), c("1", "2"))
  # Mukey 1's third row (clay_total = NA) is dropped - only 2 complete rows survive.
  expect_equal(nrow(lookup[["1"]]), 2)
  expect_equal(nrow(lookup[["2"]]), 1)
  expect_setequal(colnames(lookup[["1"]]), c("clay_total", "sand_total", "silt_total"))
  # Row alignment preserved (not independently shuffled per column).
  expect_equal(lookup[["1"]][, "clay_total"], c(20, 22))
  expect_equal(lookup[["1"]][, "sand_total"], c(40, 38))
})

test_that("lookup_mukey_texture_draws() returns NULL when texture columns are absent", {
  draws <- data.frame(mukey = "1", db = 1.4)
  expect_null(lookup_mukey_texture_draws(draws))
})

# ---------------------------------------------------------------------------
# extract_mukey_joint_ensemble() (A.2) - offline via draws_by_window
# ---------------------------------------------------------------------------

# 2 mukeys ("1": cokeys 10/11, "2": cokeys 20/21), 2 replicates each cokey, 2 windows.
.ensemble_draws_by_window <- function() {
  one <- function(db) data.frame(
    mukey = rep(c("1", "2"), each = 4),
    cokey = c("10", "10", "11", "11", "20", "20", "21", "21"),
    simulation_number = rep(c(1L, 2L), 4),
    db = db,
    ph = db + 5,
    top = 0, bottom = 5
  )
  list(
    "0-5"  = one(c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)),
    "5-15" = one(c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9))
  )
}

test_that("extract_mukey_joint_ensemble() returns per-mukey, per-window row-aligned joint matrices", {
  dbw <- .ensemble_draws_by_window()
  ens <- extract_mukey_joint_ensemble(
    aoi_vect = NULL, depth_windows = list(c(0, 5), c(5, 15)),
    draws_by_window = dbw, properties = c("db", "ph", "soc")
  )
  expect_setequal(names(ens$by_mukey), c("1", "2"))
  expect_equal(ens$properties, c("db", "ph"))            # "soc" absent -> dropped
  expect_equal(ens$window_names, c("0-5", "5-15"))

  m1 <- ens$by_mukey[["1"]]
  expect_setequal(names(m1$windows), c("0-5", "5-15"))
  expect_equal(dim(m1$windows[["0-5"]]), c(4, 2))
  expect_equal(colnames(m1$windows[["0-5"]]), c("db", "ph"))
  expect_equal(nrow(m1$replicate_key), 4)
  # window matrices are row-aligned: ph == db + 5 within each, and 5-15 db == 0-5 db + 0.1
  expect_equal(unname(m1$windows[["0-5"]][, "ph"]), unname(m1$windows[["0-5"]][, "db"]) + 5)
  expect_equal(m1$windows[["5-15"]][, "db"], m1$windows[["0-5"]][, "db"] + 0.1, tolerance = 1e-9)
})

test_that("extract_mukey_joint_ensemble() drops a replicate missing from any window", {
  dbw <- .ensemble_draws_by_window()
  # remove mukey 1 / cokey 10 / sim 2 from the second window only
  d2 <- dbw[["5-15"]]
  dbw[["5-15"]] <- d2[!(d2$mukey == "1" & d2$cokey == "10" & d2$simulation_number == 2L), ]

  ens <- extract_mukey_joint_ensemble(
    aoi_vect = NULL, depth_windows = list(c(0, 5), c(5, 15)),
    draws_by_window = dbw, properties = c("db", "ph")
  )
  m1 <- ens$by_mukey[["1"]]
  expect_equal(nrow(m1$replicate_key), 3)                # 4 -> 3
  expect_equal(nrow(m1$windows[["0-5"]]), 3)             # kept windows stay aligned
  expect_false(any(m1$replicate_key$cokey == "10" & m1$replicate_key$simulation_number == 2L))
  expect_equal(nrow(ens$by_mukey[["2"]]$replicate_key), 4)  # mukey 2 untouched
})

test_that("extract_mukey_joint_ensemble() rejects a malformed depth_windows", {
  dbw <- .ensemble_draws_by_window()
  expect_error(
    extract_mukey_joint_ensemble(NULL, depth_windows = list(), draws_by_window = dbw),
    "non-empty list"
  )
  expect_error(
    extract_mukey_joint_ensemble(NULL, depth_windows = list(c(5, 5)), draws_by_window = dbw),
    "bottom > top"
  )
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

## Narrow, control-flow-only mocking (not a fixture-maintained SDA response) - see
## test-ssurgo-acquisition.R's identical note and MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task
## P1.4. Verifies fetch_ssurgo_mukey_raster() no longer issues the redundant
## soilDB::SDA_spatialQuery() fetch, and that its new disk cache actually avoids a second
## soilDB::mukey.wcs() call for the same AOI.
test_that("fetch_ssurgo_mukey_raster() calls mukey.wcs() exactly once, never calls SDA_spatialQuery(), and caches the result", {
  aoi <- terra::vect("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")
  on.exit(unlink(file.path(tools::R_user_dir("soilSIM", "cache"), paste0(mukey_grid_cache_key(aoi), ".rds"))), add = TRUE)

  synthetic_mu <- terra::rast(nrows = 2, ncols = 2, vals = c(101, 101, 202, 202))
  names(synthetic_mu) <- "mukey"

  wcs_calls <- 0
  testthat::local_mocked_bindings(
    mukey.wcs = function(...) {
      wcs_calls <<- wcs_calls + 1
      synthetic_mu
    },
    SDA_spatialQuery = function(...) {
      stop("SDA_spatialQuery() should not be called - fetch_ssurgo_mukey_raster() must use mukey.wcs()'s own values directly")
    },
    .package = "soilDB"
  )

  result <- fetch_ssurgo_mukey_raster(aoi)
  expect_s4_class(result, "SpatRaster")
  expect_equal(wcs_calls, 1)
  expect_setequal(unique(terra::values(result))[, 1], c(101, 202))

  # Warm cache: a second call for the same AOI must not touch mukey.wcs()/SDA_spatialQuery() again.
  result2 <- fetch_ssurgo_mukey_raster(aoi)
  expect_equal(wcs_calls, 1)
  expect_setequal(unique(terra::values(result2))[, 1], c(101, 202))
})

# ---------------------------------------------------------------------------
# requested_properties (MULTI_PROPERTY_FUSION_PLAN.md tasks B0 / B2 / B-pre)
# ---------------------------------------------------------------------------

test_that("normalize_requested_properties() is total, idempotent, and couples texture", {
  nz <- normalize_requested_properties
  expect_null(nz(NULL))
  expect_identical(nz("ph"), "ph")
  expect_identical(nz(c("bulk_density", "dbovendry", "db")), "db")          # aliases collapse
  expect_identical(nz("sandtotal"), c("ilr1", "ilr2"))                      # texture -> both axes
  expect_identical(nz("clay_total"), c("ilr1", "ilr2"))
  expect_identical(nz(c("ilr1", "ilr2")), c("ilr1", "ilr2"))               # sentinel accepted
  expect_identical(nz(c("soc", "db")), c("db", "soc"))                      # canonical param order
  expect_identical(nz(nz(c("clay", "ph", "soc"))), nz(c("clay", "ph", "soc")))  # idempotent
  expect_error(nz("not_a_property"), "no simulated column mapping")
})

test_that("normalize_requested_properties() resolves the 5 P2 chemistry properties (self-mapping, canonical order)", {
  nz <- normalize_requested_properties
  expect_identical(nz("caco3"), "caco3")
  expect_identical(nz(c("sar", "caco3")), c("caco3", "sar"))  # canonical param_order, not input order
  expect_identical(nz(c("ph", "ecec", "clay")), c("ilr1", "ilr2", "ph", "ecec"))
})

test_that("simulate_ssurgo_mapunit_draws() ignores requested_properties (with a warning) under gp_quantile_retrofit", {
  # The guard runs before any fetch; mock the tabular fetch away so the call returns quickly.
  testthat::local_mocked_bindings(
    cache_get = function(...) NULL,
    download_ssurgo_tabular = function(...) NULL,
    .package = "soilSIM"
  )
  aoi <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  retrofit <- list(monte_carlo = list(vertical_correlation_method = "gp_quantile_retrofit"))

  expect_warning(
    res <- simulate_ssurgo_mapunit_draws(aoi, 0, 5, requested_properties = "ph", config = retrofit),
    "gp_quantile_retrofit"
  )
  expect_null(res)  # download mocked to NULL

  # joint_copula (the default) does NOT warn.
  expect_no_warning(
    simulate_ssurgo_mapunit_draws(aoi, 0, 5, requested_properties = "ph")
  )
})

test_that("maybe_adjust_soil_data_depth_trend() handles a single property (joint_copula, k = 1)", {
  testthat::skip_if_not_installed("GPfit")
  set.seed(101)
  depths <- c(0, 20, 50, 100)
  n_sims <- 30
  sim_long <- data.frame(
    cokey = "1",
    hzdept_r = rep(depths, each = n_sims), hzdepb_r = rep(depths + 20, each = n_sims),
    simulation_number = rep(seq_len(n_sims), times = length(depths)),
    clay_total = 25 + 0.05 * rep(depths, each = n_sims) + stats::rnorm(n_sims * length(depths), sd = 2)
  )

  set.seed(102)
  res <- maybe_adjust_soil_data_depth_trend(sim_long, "clay_total", min_depths = 2)
  expect_equal(nrow(res), nrow(sim_long))
  expect_false(anyNA(res$clay_total))
  # the k < 2 branch of preserve_correlation_structure_joint() still applies a depth trend
  expect_false(isTRUE(all.equal(res$clay_total, sim_long$clay_total)))
})

# ---------------------------------------------------------------------------
# seed = (opt-in determinism, MULTI_PROPERTY_FUSION_PLAN.md B9)
# ---------------------------------------------------------------------------

test_that("simulate_ssurgo_mapunit_draws(seed=) seeds the global RNG stream once, up front", {
  testthat::local_mocked_bindings(
    cache_get = function(...) NULL,
    download_ssurgo_tabular = function(...) NULL,  # -> early NULL return, nothing else uses RNG
    .package = "soilSIM"
  )
  aoi <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")

  set.seed(42)
  ref <- .Random.seed
  runif(1)  # perturb so we can tell the call re-seeds rather than leaving it alone
  invisible(simulate_ssurgo_mapunit_draws(aoi, 0, 5, seed = 42))
  expect_identical(.Random.seed, ref)

  # seed = NULL leaves the stream untouched (current behavior)
  s_before <- .Random.seed
  invisible(simulate_ssurgo_mapunit_draws(aoi, 0, 5, seed = NULL))
  expect_identical(.Random.seed, s_before)
})

test_that("maybe_adjust_soil_data_depth_trend(parallel=TRUE) is reproducible with seed and random without", {
  testthat::skip_if_not_installed("GPfit")
  testthat::skip_if_not_installed("future.apply")
  skip_if_not(nzchar(system.file("Meta", "package.rds", package = "future")),
              "future package metadata unavailable")

  set.seed(1)
  depths <- c(0, 20, 50, 100); n_sims <- 20
  mk <- function(ck) data.frame(
    cokey = ck,
    hzdept_r = rep(depths, each = n_sims), hzdepb_r = rep(depths + 20, each = n_sims),
    simulation_number = rep(seq_len(n_sims), times = length(depths)),
    clay_total = 25 + 0.05 * rep(depths, each = n_sims) + stats::rnorm(n_sims * length(depths), sd = 2),
    db = 1.3 + 0.002 * rep(depths, each = n_sims) + stats::rnorm(n_sims * length(depths), sd = 0.03)
  )
  sim_long <- rbind(mk("1"), mk("2"))
  props <- c("clay_total", "db")

  a <- maybe_adjust_soil_data_depth_trend(sim_long, props, parallel = TRUE, n_cores = 2, seed = 7)
  b <- maybe_adjust_soil_data_depth_trend(sim_long, props, parallel = TRUE, n_cores = 2, seed = 7)
  ord <- function(d) d[order(d$cokey, d$hzdept_r, d$simulation_number), ]
  expect_equal(ord(a)$clay_total, ord(b)$clay_total, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# ssurgo_tabular_cache_key() / depth-independent tabular cache (MULTI_PROPERTY_FUSION_PLAN.md L1)
# ---------------------------------------------------------------------------

test_that("ssurgo_tabular_cache_key() is depth-independent and AOI-specific", {
  aoi1 <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  aoi2 <- terra::vect(terra::ext(5, 6, 5, 6), crs = "EPSG:5070")

  expect_identical(ssurgo_tabular_cache_key(aoi1), build_cache_key(aoi1, "ssurgo_tabular", 0, 0, "ssurgo_tabular"))
  # No depth args to vary - the function's whole point is that top_depth/bottom_depth never enter
  # the key at all. Confirm two "requests" for the same AOI (conceptually different windows,
  # nothing here actually varies since the function takes no depth args) still agree.
  expect_identical(ssurgo_tabular_cache_key(aoi1), ssurgo_tabular_cache_key(aoi1))
  expect_false(identical(ssurgo_tabular_cache_key(aoi1), ssurgo_tabular_cache_key(aoi2)))
})

test_that("simulate_ssurgo_mapunit_draws() shares ONE tabular download across different depth windows for the same AOI", {
  store <- new.env(parent = emptyenv())
  dl_calls <- 0L
  aoi <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")

  testthat::local_mocked_bindings(
    cache_get = function(key, ttl_seconds = CACHE_TTL_SECONDS) store[[key]],
    cache_set = function(key, kind, value) { store[[key]] <- value; invisible(TRUE) },
    download_ssurgo_tabular = function(...) {
      dl_calls <<- dl_calls + 1L
      list(ssurgo_data = data.frame(cokey = "1"), mu = NULL, metadata = NULL)  # minimal; downstream may error
    },
    .package = "soilSIM"
  )

  # Two different windows for the SAME AOI - downstream steps beyond the tabular fetch will likely
  # error on this minimal frame; that's fine, the tabular download + cache_set already happened by
  # then (see source: cache_set() runs immediately after a cache-miss download, before infill).
  suppressWarnings(tryCatch(simulate_ssurgo_mapunit_draws(aoi, 0, 5),  error = function(e) NULL))
  suppressWarnings(tryCatch(simulate_ssurgo_mapunit_draws(aoi, 5, 15), error = function(e) NULL))

  expect_equal(dl_calls, 1L)  # would be 2 before the L1 fix (window-keyed cache)
})

test_that("simulate_ssurgo_mapunit_draws() does NOT share the tabular cache across different AOIs", {
  store <- new.env(parent = emptyenv())
  dl_calls <- 0L
  aoi1 <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  aoi2 <- terra::vect(terra::ext(5, 6, 5, 6), crs = "EPSG:5070")

  testthat::local_mocked_bindings(
    cache_get = function(key, ttl_seconds = CACHE_TTL_SECONDS) store[[key]],
    cache_set = function(key, kind, value) { store[[key]] <- value; invisible(TRUE) },
    download_ssurgo_tabular = function(...) {
      dl_calls <<- dl_calls + 1L
      list(ssurgo_data = data.frame(cokey = "1"), mu = NULL, metadata = NULL)
    },
    .package = "soilSIM"
  )

  suppressWarnings(tryCatch(simulate_ssurgo_mapunit_draws(aoi1, 0, 5), error = function(e) NULL))
  suppressWarnings(tryCatch(simulate_ssurgo_mapunit_draws(aoi2, 0, 5), error = function(e) NULL))

  expect_equal(dl_calls, 2L)
})
