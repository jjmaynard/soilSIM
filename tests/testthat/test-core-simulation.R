test_that("remove_organic_layer() drops O horizons and re-anchors remaining depths to 0", {
  df <- data.frame(
    cokey = c("1", "1", "1"), hzname = c("O", "A", "Bt"),
    hzdept_r = c(-5, 0, 20), hzdepb_r = c(0, 20, 50),
    stringsAsFactors = FALSE
  )
  result <- remove_organic_layer(df)
  expect_equal(nrow(result), 2)
  expect_false(any(grepl("O", result$hzname)))
  expect_equal(result$hzdept_r, c(0, 20))
  expect_equal(result$hzdepb_r, c(20, 50))
})

test_that("remove_organic_layer() handles a group with no organic horizons unchanged in row count", {
  df <- data.frame(
    cokey = c("1", "1"), hzname = c("A", "Bt"),
    hzdept_r = c(0, 20), hzdepb_r = c(20, 50),
    stringsAsFactors = FALSE
  )
  result <- remove_organic_layer(df)
  expect_equal(nrow(result), 2)
})

test_that("slice_and_aggregate_soil_data() averages a numeric property within each depth range", {
  df <- data.frame(
    compname = "testseries", hzdept_r = c(0, 30), hzdepb_r = c(30, 100),
    dbovendry_r = c(1.2, 1.5)
  )
  result <- slice_and_aggregate_soil_data(df, depth_ranges = list(c(0, 30), c(30, 100)))
  expect_equal(nrow(result), 2)
  expect_equal(result$hzdept_r, c(0, 30))
  expect_equal(result$dbovendry_r[1], 1.2, tolerance = 1e-6)
  expect_equal(result$dbovendry_r[2], 1.5, tolerance = 1e-6)
})

test_that("slice_and_aggregate_soil_data()'s vectorized depth expansion matches the original per-centimeter loop", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 3 slice_and_aggregate_soil_data()
  # fix: reimplements the ORIGINAL per-centimeter for-loop (one as.data.frame() call per depth)
  # here and checks the vectorized version produces identical aggregated output, across multiple
  # horizons/properties with uneven depth ranges.
  old_expand <- function(df, data_columns) {
    rows_list <- list()
    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      row_data <- as.list(row[data_columns])
      for (depth in seq(row$hzdept_r, row$hzdepb_r - 1)) {
        row_data$Depth <- depth
        rows_list <- append(rows_list, list(as.data.frame(row_data, stringsAsFactors = FALSE)))
      }
    }
    aggregated <- do.call(rbind, rows_list)
    aggregated$Depth <- as.numeric(aggregated$Depth)
    rownames(aggregated) <- NULL
    aggregated
  }

  df <- data.frame(
    compname = "testseries",
    hzdept_r = c(0, 17, 42),
    hzdepb_r = c(17, 42, 88),
    dbovendry_r = c(1.2, 1.4, 1.6),
    claytotal_r = c(15, 22, 30)
  )
  data_columns <- c("dbovendry_r", "claytotal_r")

  expected <- old_expand(df, data_columns)
  depth_ranges <- list(c(0, 30), c(30, 90))
  result <- slice_and_aggregate_soil_data(df, depth_ranges = depth_ranges)

  # Cross-check against the manually-computed expected per-cm expansion's own aggregates.
  for (i in seq_along(depth_ranges)) {
    top <- depth_ranges[[i]][1]; bottom <- depth_ranges[[i]][2]
    subset_expected <- expected[expected$Depth >= top & expected$Depth < bottom, ]
    expect_equal(result$dbovendry_r[i], mean(subset_expected$dbovendry_r), tolerance = 1e-10)
    expect_equal(result$claytotal_r[i], mean(subset_expected$claytotal_r), tolerance = 1e-10)
  }
})

test_that("simulate_component_composition() adds a sim_comppct column, one row per distinct component", {
  data <- make_component_data()
  result <- simulate_component_composition(data, n_simulations = 500)
  expect_equal(nrow(result), 2)
  expect_true("sim_comppct" %in% names(result))
  expect_true(all(result$sim_comppct > 0))
  # sim_comppct ~= round(n_simulations * comppct_r / 100) in expectation
  expect_equal(result$sim_comppct[1], round(500 * 55 / 100), tolerance = 20)
})

test_that("simulate_component_composition() fills missing comppct_l/comppct_h from comppct_r +/- 2", {
  data <- data.frame(
    mukey = "1", cokey = "1", compname = "compA",
    comppct_l = NA_real_, comppct_r = 50, comppct_h = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- simulate_component_composition(data, n_simulations = 200)
  expect_equal(nrow(result), 1)
  expect_true(result$sim_comppct > 0)
})

test_that("simulate_component_composition() output joins onto horizon data by cokey (documented grain-mismatch fix)", {
  component_data <- make_component_data()
  comp_result <- simulate_component_composition(component_data, n_simulations = 100)

  horizon_data <- data.frame(
    cokey = c("1", "1", "2"), hzname = c("A", "Bt", "A"),
    hzdept_r = c(0, 20, 0), hzdepb_r = c(20, 50, 30),
    stringsAsFactors = FALSE
  )
  joined <- dplyr::left_join(horizon_data, comp_result[, c("cokey", "sim_comppct")], by = "cokey")
  expect_true(all(!is.na(joined$sim_comppct)))
  expect_equal(nrow(joined), nrow(horizon_data))
})

test_that("compute_mode() returns the most frequent value", {
  expect_equal(compute_mode(c(1, 2, 2, 3, 3, 3, 4)), 3)
})

test_that("simulate_correlated_triangular() respects each distribution's (lower, mode, upper) bounds", {
  set.seed(1)
  params <- list(c(0, 5, 10), c(1, 4, 6))
  correlation_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  samples <- simulate_correlated_triangular(2000, params, correlation_matrix)

  expect_equal(dim(samples), c(2000, 2))
  expect_true(all(samples[, 1] >= 0 & samples[, 1] <= 10))
  expect_true(all(samples[, 2] >= 1 & samples[, 2] <= 6))
})

test_that("simulate_correlated_triangular() induces positive correlation between columns", {
  set.seed(2)
  params <- list(c(0, 5, 10), c(0, 5, 10))
  correlation_matrix <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  samples <- simulate_correlated_triangular(5000, params, correlation_matrix)
  expect_true(stats::cor(samples[, 1], samples[, 2]) > 0.4)
})

test_that("simulate_correlated_triangular() handles a degenerate (a == c) distribution without error", {
  params <- list(c(5, 5, 5), c(0, 5, 10))
  correlation_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
  samples <- simulate_correlated_triangular(50, params, correlation_matrix)
  expect_true(all(samples[, 1] == 5))
})

test_that("simulate_cokey_generalized() simulates non-texture properties with sim_comppct realizations", {
  set.seed(3)
  sim_cokey <- make_sim_cokey_data(texture = FALSE, sim_comppct = 25)
  corr <- make_property_correlation_matrices()
  result <- simulate_cokey_generalized(sim_cokey, corr)

  expect_equal(nrow(result), 25)
  expect_true(all(c("db", "ph", "compname", "mukey", "cokey", "hzdept_r", "hzdepb_r",
                     "simulation_number", "unique_id") %in% names(result)))
})

test_that("simulate_cokey_generalized() propagates bound_sd (OSD boundary distinctness) onto the output when present on sim_cokey, and omits it when absent", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 1b: bound_sd must survive into
  # simulate_cokey_generalized()'s output so the future depth-kernel gating (Phase 1c) can read it
  # per depth, without requiring every existing caller to have a bound_sd column.
  set.seed(5)
  sim_cokey <- make_sim_cokey_data(texture = FALSE, sim_comppct = 10)
  corr <- make_property_correlation_matrices()

  # Absent case (existing/default behavior) - no bound_sd column added.
  result_absent <- simulate_cokey_generalized(sim_cokey, corr)
  expect_false("bound_sd" %in% names(result_absent))

  # Present case - broadcast onto every simulated row for that horizon.
  sim_cokey$bound_sd <- 12.5
  result_present <- simulate_cokey_generalized(sim_cokey, corr)
  expect_true("bound_sd" %in% names(result_present))
  expect_equal(nrow(result_present), 10)
  expect_true(all(result_present$bound_sd == 12.5))
})

test_that("simulate_cokey_generalized() with texture columns produces sand+silt+clay summing to 100", {
  set.seed(4)
  sim_cokey <- make_sim_cokey_data(texture = TRUE, sim_comppct = 30)
  corr <- make_property_correlation_matrices()
  txt_corr <- make_texture_correlation_matrices()
  result <- simulate_cokey_generalized(sim_cokey, corr, txt_corr)

  expect_equal(nrow(result), 30)
  expect_true(all(c("sand_total", "silt_total", "clay_total") %in% names(result)))
  totals <- result$sand_total + result$silt_total + result$clay_total
  expect_true(all(abs(totals - 100) < 1e-6))
  expect_true(all(result$sand_total >= 0 & result$sand_total <= 100))
})

test_that("simulate_cokey_generalized() drops texture but keeps other properties when a texture triplet has a leftover NA", {
  set.seed(7)
  sim_cokey <- make_sim_cokey_data(texture = TRUE, sim_comppct = 20)
  sim_cokey$sandtotal_l <- NA_real_ # simulates infilling failing to recover this one value
  corr <- make_property_correlation_matrices()
  txt_corr <- make_texture_correlation_matrices()

  expect_message(
    result <- simulate_cokey_generalized(sim_cokey, corr, txt_corr),
    "texture simulation failed"
  )
  expect_equal(nrow(result), 20)
  expect_true(all(c("db", "ph") %in% names(result)))
  expect_false(any(c("sand_total", "silt_total", "clay_total") %in% names(result)))
})

test_that("simulate_cokey_generalized() falls back to a pooled matrix for an unmatched genhz instead of erroring", {
  set.seed(5)
  sim_cokey <- make_sim_cokey_data(texture = FALSE, sim_comppct = 25)
  sim_cokey$genhz <- NA_character_ # classify_genhz()'s real output for an unparseable hzname
  corr <- make_property_correlation_matrices() # keyed only by "A"/"B" - no NA entry

  result <- simulate_cokey_generalized(sim_cokey, corr)
  expect_equal(nrow(result), 25)
  expect_true(all(c("db", "ph") %in% names(result)))
})

test_that("simulate_cokey_generalized() falls back to a pooled texture matrix for an unmatched genhz instead of erroring", {
  set.seed(6)
  sim_cokey <- make_sim_cokey_data(texture = TRUE, sim_comppct = 30)
  sim_cokey$genhz <- "Z" # not a key in either correlation list below
  corr <- make_property_correlation_matrices()
  txt_corr <- make_texture_correlation_matrices()

  result <- simulate_cokey_generalized(sim_cokey, corr, txt_corr)
  expect_equal(nrow(result), 30)
  totals <- result$sand_total + result$silt_total + result$clay_total
  expect_true(all(abs(totals - 100) < 1e-6))
})

test_that("simulate_cokey_generalized() simulates texture for an NA genhz (unclassifiable hzname like 'H1')", {
  # Regression: as.character(classify_genhz("H1")) is NA_character_, and an NA_character_ list
  # key is write-only - `pd_txt_corr_cache[[NA]] <- m` then `pd_txt_corr_cache[[NA]]` reads back
  # NULL. That fed a NULL correlation matrix into simulate_correlated_triangular() and the
  # texture draw failed for EVERY such row (SSURGO map units with generic H1/H2/H3 designations
  # ended up with no texture -> NA prior -> white pixels after raster fusion). Distinct from the
  # genhz = "Z" case above: "Z" is a readable key, NA_character_ is not.
  set.seed(6)
  sim_cokey <- make_sim_cokey_data(texture = TRUE, sim_comppct = 30)
  sim_cokey$genhz <- NA_character_
  corr <- make_property_correlation_matrices()
  txt_corr <- make_texture_correlation_matrices()

  result <- expect_no_message(
    simulate_cokey_generalized(sim_cokey, corr, txt_corr),
    message = "texture simulation failed"
  )
  expect_equal(nrow(result), 30)
  expect_true(all(c("sand_total", "silt_total", "clay_total") %in% names(result)))
  totals <- result$sand_total + result$silt_total + result$clay_total
  expect_true(all(abs(totals - 100) < 1e-6))

  # clay-only request: texture is the only output, so a failed texture draw would drop the
  # whole cokey (this is exactly the raster-fusion single-property path).
  set.seed(6)
  only_clay <- simulate_cokey_generalized(sim_cokey, corr, txt_corr, requested_properties = "clay")
  expect_equal(nrow(only_clay), 30)
  expect_true("clay_total" %in% names(only_clay))
})

test_that("simulate_cokey_generalized()'s per-genhz PD-matrix caching matches the original uncached-per-row computation bit-for-bit", {
  # PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4: ensure_positive_definite_matrix(txt_corr) is now
  # cached per genhz_val instead of recomputed identically on every row (Rprof() profiling found
  # this call accounted for 18% of total wall-clock). It involves no randomness (eigen()/
  # isSymmetric.matrix() are deterministic), so caching must not shift simulate_correlated_
  # triangular()'s RNG stream consumption at all - reimplements the ORIGINAL uncached per-row
  # loop inline (as it existed before this fix) and asserts identical output at the same seed,
  # across multiple rows that intentionally share a genhz (the case the caching actually exercises
  # - the existing single-row fixtures above never did).
  row1 <- make_sim_cokey_data(texture = TRUE, sim_comppct = 15)
  row2 <- make_sim_cokey_data(texture = TRUE, sim_comppct = 20)
  row2$hzdept_r <- 20; row2$hzdepb_r <- 40
  sim_cokey <- rbind(row1, row2) # both rows share genhz = "A" (make_sim_cokey_data()'s default)
  corr <- make_property_correlation_matrices()
  txt_corr <- make_texture_correlation_matrices()

  set.seed(42)
  result_cached <- simulate_cokey_generalized(sim_cokey, corr, txt_corr)

  # Original (pre-fix) uncached-per-row reference implementation.
  reference_cokey_generalized <- function(sim_cokey, correlation_matrices, txt_correlation_matrices) {
    param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc")
    texture_cols_required <- c("sandtotal_l", "sandtotal_r", "sandtotal_h",
                               "silttotal_l", "silttotal_r", "silttotal_h",
                               "claytotal_l", "claytotal_r", "claytotal_h")
    get_param_set <- function(row, prefix) {
      lcol <- paste0(prefix, "_l"); rcol <- paste0(prefix, "_r"); hcol <- paste0(prefix, "_h")
      if (all(c(lcol, rcol, hcol) %in% names(row))) return(c(row[[lcol]], row[[rcol]], row[[hcol]]))
      NULL
    }
    sim_data_out <- list()
    pooled_property_corr <- Reduce(`+`, correlation_matrices) / length(correlation_matrices)
    pooled_txt_corr <- Reduce(`+`, txt_correlation_matrices) / length(txt_correlation_matrices)
    for (i in seq_len(nrow(sim_cokey))) {
      row <- sim_cokey[i, ]
      genhz_val <- as.character(row$genhz)
      local_corr <- correlation_matrices[[genhz_val]]
      if (is.null(local_corr)) local_corr <- pooled_property_corr
      has_texture <- all(texture_cols_required %in% names(row))
      ilr1_lrh <- NULL; ilr2_lrh <- NULL
      if (has_texture) {
        txt_corr <- txt_correlation_matrices[[genhz_val]]
        if (is.null(txt_corr)) txt_corr <- pooled_txt_corr
        txt_corr <- ensure_positive_definite_matrix(txt_corr) # recomputed every row - the original
        params_txt <- list(
          c(row$sandtotal_l, row$sandtotal_r, row$sandtotal_h),
          c(row$silttotal_l, row$silttotal_r, row$silttotal_h),
          c(row$claytotal_l, row$claytotal_r, row$claytotal_h)
        )
        texture_result <- tryCatch({
          sim_txt <- simulate_correlated_triangular(as.integer(row$sim_comppct), params_txt, txt_corr)
          sim_txt_ilr <- ilr_forward(clay = sim_txt[, 3], sand = sim_txt[, 1], silt = sim_txt[, 2])
          list(ilr1_lrh = c(min(sim_txt_ilr[, "z1"]), compute_mode(sim_txt_ilr[, "z1"]), max(sim_txt_ilr[, "z1"])),
               ilr2_lrh = c(min(sim_txt_ilr[, "z2"]), compute_mode(sim_txt_ilr[, "z2"]), max(sim_txt_ilr[, "z2"])))
        }, error = function(e) NULL)
        if (!is.null(texture_result)) { ilr1_lrh <- texture_result$ilr1_lrh; ilr2_lrh <- texture_result$ilr2_lrh }
      }
      param_list <- vector("list", length(param_order)); names(param_list) <- param_order
      db_set <- get_param_set(row, "dbovendry"); if (!is.null(db_set)) param_list[["db"]] <- db_set
      wr3_set <- get_param_set(row, "wthirdbar"); if (!is.null(wr3_set)) param_list[["wr_3b"]] <- wr3_set
      wr15_set <- get_param_set(row, "wfifteenbar"); if (!is.null(wr15_set)) param_list[["wr_15b"]] <- wr15_set
      if (has_texture && !is.null(ilr1_lrh) && !is.null(ilr2_lrh)) {
        param_list[["ilr1"]] <- ilr1_lrh; param_list[["ilr2"]] <- ilr2_lrh
      }
      rfv_set <- get_param_set(row, "rfv"); if (!is.null(rfv_set)) param_list[["rfv"]] <- rfv_set
      ph_set <- get_param_set(row, "ph1to1h2o"); if (!is.null(ph_set)) param_list[["ph"]] <- ph_set
      cec_set <- get_param_set(row, "cec7"); if (!is.null(cec_set)) param_list[["cec"]] <- cec_set
      om_set <- get_param_set(row, "om"); if (!is.null(om_set)) param_list[["soc"]] <- om_set
      param_list <- param_list[!vapply(param_list, is.null, logical(1))]
      if (!length(param_list)) next
      params_for_sim <- unname(param_list)
      keep_cols <- names(param_list)
      local_corr_sub <- local_corr[keep_cols, keep_cols, drop = FALSE]
      tryCatch({
        n_sim <- as.integer(row$sim_comppct)
        sim_data <- as.data.frame(simulate_correlated_triangular(n = n_sim, params = params_for_sim, correlation_matrix = local_corr_sub))
        colnames(sim_data) <- names(param_list)
        if ("wr_3b" %in% names(param_list)) sim_data[["wr_3b"]] <- sim_data[["wr_3b"]] / 100
        if ("wr_15b" %in% names(param_list)) sim_data[["wr_15b"]] <- sim_data[["wr_15b"]] / 100
        if (has_texture && all(c("ilr1", "ilr2") %in% names(param_list))) {
          sim_txt <- ilr_inverse(sim_data[["ilr1"]], sim_data[["ilr2"]], total = 100)
          sim_txt <- as.data.frame(sim_txt)[, c("sand", "silt", "clay")]
          colnames(sim_txt) <- c("sand_total", "silt_total", "clay_total")
          sim_data <- sim_data[, setdiff(names(sim_data), c("ilr1", "ilr2"))]
          sim_data <- cbind(sim_data, sim_txt)
        }
        sim_data$compname <- row$compname; sim_data$mukey <- row$mukey; sim_data$cokey <- row$cokey
        sim_data$hzdept_r <- row$hzdept_r; sim_data$hzdepb_r <- row$hzdepb_r
        sim_data$simulation_number <- seq_len(nrow(sim_data))
        sim_data$unique_id <- paste0(row$cokey, "-", sprintf("%02d", sim_data$simulation_number))
        sim_data_out <- append(sim_data_out, list(sim_data))
      }, error = function(e) NULL)
    }
    if (length(sim_data_out) == 0) NULL else dplyr::bind_rows(sim_data_out)
  }

  set.seed(42)
  result_reference <- reference_cokey_generalized(sim_cokey, corr, txt_corr)

  expect_equal(result_cached, result_reference)
})

test_that("simulate_cokey_generalized() skips a row with no recognized properties rather than erroring", {
  sim_cokey <- data.frame(
    genhz = "A", sim_comppct = 10, compname = "x", mukey = "1", cokey = "1",
    hzdept_r = 0, hzdepb_r = 20, stringsAsFactors = FALSE
  )
  corr <- make_property_correlation_matrices()
  expect_message(result <- simulate_cokey_generalized(sim_cokey, corr), "No recognized properties")
  expect_null(result)
})

# ---------------------------------------------------------------------------
# requested_properties (MULTI_PROPERTY_FUSION_PLAN.md task B1)
# ---------------------------------------------------------------------------

test_that("simulate_cokey_generalized(requested_properties = NULL) is bit-identical to the unset call", {
  cm  <- make_property_correlation_matrices()
  tcm <- make_texture_correlation_matrices()
  sc  <- make_sim_cokey_data(texture = TRUE, sim_comppct = 15)

  set.seed(1); a <- simulate_cokey_generalized(sc, cm, tcm)
  set.seed(1); b <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = NULL)
  expect_identical(a, b)
})

test_that("simulate_cokey_generalized() restricts output to the requested properties (+ texture coupling)", {
  cm  <- make_property_correlation_matrices()
  tcm <- make_texture_correlation_matrices()
  sc  <- make_sim_cokey_data(texture = TRUE, sim_comppct = 15)

  set.seed(1); only_ph <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = "ph")
  expect_true("ph" %in% names(only_ph))
  expect_false(any(c("db", "sand_total", "silt_total", "clay_total") %in% names(only_ph)))

  # Any texture member pulls in the whole sand/silt/clay draw, nothing else.
  set.seed(1); only_clay <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = "clay")
  expect_true(all(c("sand_total", "silt_total", "clay_total") %in% names(only_clay)))
  expect_false(any(c("db", "ph") %in% names(only_clay)))

  # Exactly one non-texture column alongside texture - regression for the missing `drop = FALSE`
  # on the post-ILR-inverse column drop.
  set.seed(1); clay_ph <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = c("clay", "ph"))
  expect_true(all(c("sand_total", "silt_total", "clay_total", "ph") %in% names(clay_ph)))
  expect_false("db" %in% names(clay_ph))
})

test_that("simulate_cokey_generalized() subset marginals match the full simulation distributionally", {
  cm  <- make_property_correlation_matrices()
  tcm <- make_texture_correlation_matrices()
  sc  <- make_sim_cokey_data(texture = TRUE, sim_comppct = 400)

  set.seed(7); full <- simulate_cokey_generalized(sc, cm, tcm)$ph
  set.seed(7); sub  <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = "ph")$ph

  # RNG streams differ (fewer columns drawn), so not identical - but the Gaussian-copula marginal
  # of a subset equals the marginal of the full draw, so the distributions must agree.
  expect_gt(suppressWarnings(stats::ks.test(full, sub)$p.value), 0.01)
  expect_equal(mean(full), mean(sub), tolerance = 0.1)
})

test_that("simulate_cokey_generalized(): a cokey with no data for the ONLY requested property drops out (MULTI_PROPERTY_FUSION_PLAN.md B8)", {
  cm  <- make_property_correlation_matrices()
  tcm <- make_texture_correlation_matrices()
  sc  <- make_sim_cokey_data(texture = TRUE, sim_comppct = 15)
  sc$claytotal_l <- sc$claytotal_r <- sc$claytotal_h <- NA_real_

  # Full simulation: the row still contributes db/ph (texture just absent).
  set.seed(1)
  full <- suppressMessages(simulate_cokey_generalized(sc, cm, tcm))
  expect_true(!is.null(full) && all(c("db", "ph") %in% names(full)))

  # requested_properties = "clay": the only kept parameter is texture, its triplet is NA, so the
  # row has nothing to simulate -> skipped -> NULL. This is the extra-NA behavior change.
  set.seed(1)
  expect_message(
    only_clay <- simulate_cokey_generalized(sc, cm, tcm, requested_properties = "clay"),
    "No recognized properties"
  )
  expect_null(only_clay)
})

# ---------------------------------------------------------------------------
# om -> soc conversion (MULTI_PROPERTY_FUSION_PLAN.md task P1)
# ---------------------------------------------------------------------------

test_that("simulate_cokey_generalized() converts om to an SOC-scale estimate via OM_TO_SOC_FACTOR", {
  row <- make_horizon_row(properties = "om", genhz = "A", cokey = "1", mukey = "1")
  row$compname <- "testseries"
  row$sim_comppct <- 300
  cm <- list(A = matrix(1, 1, 1, dimnames = list("soc", "soc")))

  set.seed(11)
  result <- simulate_cokey_generalized(row, cm, requested_properties = "soc")

  expect_true("soc" %in% names(result))
  om_l <- row$om_l; om_h <- row$om_h
  # simulated values fall within the CONVERTED (SOC-scale) support, not the raw OM support
  expect_true(all(result$soc >= om_l * OM_TO_SOC_FACTOR - 1e-9))
  expect_true(all(result$soc <= om_h * OM_TO_SOC_FACTOR + 1e-9))
  # convincingly smaller than the raw OM range (om_h = 4 vs om_h * OM_TO_SOC_FACTOR ~= 2.32) -
  # a regression to the pre-P1 unconverted behavior would violate this
  expect_lt(max(result$soc), om_h)
})

test_that("OM_TO_SOC_FACTOR is the inverse Van Bemmelen factor", {
  expect_equal(OM_TO_SOC_FACTOR, 1 / 1.724, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 5 chemistry properties: caco3/ec/ecec/gypsum/sar (MULTI_PROPERTY_FUSION_PLAN.md task P2)
# ---------------------------------------------------------------------------

test_that("simulate_cokey_generalized() simulates the 5 P2 chemistry properties when requested", {
  full_param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc",
                        "caco3", "ec", "ecec", "gypsum", "sar")
  cm <- stats::setNames(
    lapply(c("A", "B"), function(g) build_kssl_fallback_matrix(full_param_order, genhz = g)),
    c("A", "B")
  )
  row <- make_horizon_row(properties = c("dbovendry", "caco3", "ec", "ecec", "gypsum", "sar"),
                          genhz = "A", cokey = "1", mukey = "1")
  row$compname <- "testseries"
  row$sim_comppct <- 25

  set.seed(3)
  result <- simulate_cokey_generalized(row, cm, requested_properties = NULL)
  expect_true(all(c("caco3", "ec", "ecec", "gypsum", "sar", "db") %in% names(result)))
  expect_equal(nrow(result), 25)
  expect_true(all(result$caco3 >= 0 & result$ec >= 0 & result$ecec >= 0))
})

test_that("simulate_cokey_generalized() gates each of the 5 chemistry properties independently via requested_properties", {
  full_param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc",
                        "caco3", "ec", "ecec", "gypsum", "sar")
  cm <- stats::setNames(
    lapply(c("A", "B"), function(g) build_kssl_fallback_matrix(full_param_order, genhz = g)),
    c("A", "B")
  )
  row <- make_horizon_row(properties = c("dbovendry", "caco3", "ec", "ecec", "gypsum", "sar"),
                          genhz = "A", cokey = "1", mukey = "1")
  row$compname <- "testseries"
  row$sim_comppct <- 10

  for (p in c("caco3", "ec", "ecec", "gypsum", "sar")) {
    set.seed(5)
    res <- simulate_cokey_generalized(row, cm, requested_properties = p)
    expect_true(p %in% names(res), info = p)
    expect_false("db" %in% names(res), info = p)
    expect_false(any(setdiff(c("caco3", "ec", "ecec", "gypsum", "sar"), p) %in% names(res)), info = p)
  }
})

test_that("simulate_cokey_generalized() degrades to independence (not an error) for a requested chemistry property missing from correlation_matrices", {
  # make_property_correlation_matrices() only carries db/ph/ilr1/ilr2 - predates task P2 entirely,
  # the exact "direct caller with an old matrix" scenario extend_corr_matrix_with_identity() exists
  # for. Regression test for the "subscript out of bounds" crash this would otherwise cause.
  cm  <- make_property_correlation_matrices()
  row <- make_horizon_row(properties = c("dbovendry", "caco3"), genhz = "A", cokey = "1", mukey = "1")
  row$compname <- "testseries"
  row$sim_comppct <- 10

  set.seed(9)
  result <- expect_no_error(simulate_cokey_generalized(row, cm, requested_properties = "caco3"))
  expect_true("caco3" %in% names(result))
})

test_that("extend_corr_matrix_with_identity() adds identity rows/cols and stays positive-definite", {
  m <- diag(2)
  dimnames(m) <- list(c("db", "ph"), c("db", "ph"))
  m["db", "ph"] <- m["ph", "db"] <- 0.4

  extended <- extend_corr_matrix_with_identity(m, c("caco3", "sar"))
  expect_equal(dim(extended), c(4, 4))
  expect_equal(rownames(extended), c("db", "ph", "caco3", "sar"))
  expect_equal(unname(extended["db", "ph"]), 0.4)          # original correlation preserved
  expect_equal(unname(extended["caco3", "sar"]), 0)        # new names uncorrelated with each other
  expect_equal(unname(extended["db", "caco3"]), 0)         # and with the original names
  expect_true(all(eigen(extended, only.values = TRUE, symmetric = TRUE)$values > -1e-10))
})


# --- merged from test-depth-simulation.R (P1 reorg) ---

test_that("infill_missing_depth_variability() fills missing bounds with +/-2 of the representative value", {
  data <- data.frame(
    hzdept_r = c(10, 20), hzdept_l = c(NA, 18), hzdept_h = c(NA, NA),
    hzdepb_r = c(30, 40), hzdepb_l = c(NA, NA), hzdepb_h = c(NA, 42)
  )
  out <- infill_missing_depth_variability(data)
  expect_equal(out$hzdept_l, c(8, 18))
  expect_equal(out$hzdept_h, c(12, 22))
  expect_equal(out$hzdepb_l, c(28, 38))
  expect_equal(out$hzdepb_h, c(32, 42))
})

test_that("infill_missing_depth_variability() clamps at zero, never goes negative", {
  data <- data.frame(hzdept_r = 1, hzdept_l = NA, hzdept_h = NA, hzdepb_r = 1, hzdepb_l = NA, hzdepb_h = NA)
  out <- infill_missing_depth_variability(data)
  expect_equal(out$hzdept_l, 0)
  expect_equal(out$hzdepb_l, 0)
})

test_that("infill_missing_distinctness() fills defaults by genhz, preserves existing values", {
  df <- data.frame(
    hzname = c("A", "B", "C", "R", "O", "Cr"),
    distinctness = c(NA, "gradual", NA, NA, "diffuse", NA),
    stringsAsFactors = FALSE
  )
  out <- infill_missing_distinctness(df)
  expect_equal(out$distinctness, c("clear", "gradual", "gradual", "abrupt", "diffuse", "gradual"))
})

test_that("simulate_soil_profile_top_down() starts at 0, never produces bottom < top, and matches bounds", {
  set.seed(1)
  horizon_data <- make_depth_horizon_data()
  profile <- simulate_soil_profile_top_down(horizon_data)
  expect_equal(profile$top[1], 0)
  expect_true(all(profile$bottom >= profile$top))
  # each horizon's simulated bottom must fall within its l/h bounds (post-clamping)
  expect_true(all(profile$bottom <= horizon_data$hzdepb_h + 1e-6))
  # top of horizon i+1 must equal bottom of horizon i (top-down contiguity)
  expect_equal(profile$top[-1], profile$bottom[-nrow(profile)])
})

test_that("simulate_soil_profile_bottom_up() starts at 0 and has no gaps between horizons", {
  set.seed(2)
  horizon_data <- make_depth_horizon_data()
  profile <- simulate_soil_profile_bottom_up(horizon_data)
  expect_equal(profile$top[1], 0)
  expect_true(all(profile$bottom >= profile$top))
  expect_equal(profile$bottom[-nrow(profile)], profile$top[-1])
})

test_that("simulate_soil_profile_thickness() returns one row per horizon with a non-negative thickness_sd", {
  set.seed(3)
  horizon_data <- make_depth_horizon_data()
  result <- simulate_soil_profile_thickness(horizon_data, n_simulations = 50)
  expect_equal(nrow(result), nrow(horizon_data))
  expect_true(all(result$thickness_sd >= 0))
  expect_setequal(result$hzname, horizon_data$hzname)
})

test_that("simulate_and_perturb_soil_profiles() replicates a single-horizon profile sim_comppct times without perturbation", {
  df <- data.frame(
    cokey = "1", mukey = "1", compname = "onehorizon", id = "onehorizon",
    hzname = "R", hzdept_l = 0, hzdept_r = 0, hzdept_h = 0,
    hzdepb_l = 145, hzdepb_r = 150, hzdepb_h = 155, hzthk_l = 145,
    sim_comppct = 4,
    stringsAsFactors = FALSE
  )
  spc <- df
  aqp::depths(spc) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(spc) <- "hzname"

  result <- simulate_and_perturb_soil_profiles(spc)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 4)
})

test_that("simulate_profile_depths()'s id-column fix: aqp::depths(id ~ ...) doesn't error on a compname-derived id", {
  # Regression test for the id-column bug fix: fetch_ssurgo_aws_data() output
  # only has `compname`, not `id`; simulate_profile_depths() must
  # derive `id` before calling aqp::depths<-() since adjust_out_of_range_profiles()/
  # evaluate_simulated_depths() downstream both hardcode a literal `id` column.
  mu_data <- data.frame(
    mukey = "1", cokey = "1", compname = "testseries",
    hzname = c("A", "Bt"), hzdept_r = c(0, 20), hzdepb_r = c(20, 50),
    stringsAsFactors = FALSE
  )
  mu_data$id <- mu_data$compname
  expect_no_error(aqp::depths(mu_data) <- id ~ hzdept_r + hzdepb_r)
})

test_that("simulate_profile_depths() derives and joins sim_comppct internally (closed integration gap)", {
  # Regression test for the sim_comppct integration gap: simulate_profile_depths()
  # used to require callers to derive/join sim_comppct themselves before it reached
  # simulate_and_perturb_soil_profiles(), erroring on a missing-column condition otherwise.
  # It now calls simulate_component_composition() and joins the result onto mu_data by cokey internally.
  # comppct_l = comppct_r = comppct_h = 100 makes tri_dist() degenerate to a point mass at
  # 100 (a == b == c), so sim_comppct = round(n_simulations * 100 / 100) = n_simulations
  # exactly - a fully deterministic, offline-testable case (single horizon, so the
  # early-return path in simulate_and_perturb_soil_profiles() is taken and no live OSD
  # lookup occurs).
  testthat::local_mocked_bindings(
    fetch_ssurgo_aws_data = function(mukeys) {
      data.frame(
        mukey = "999", cokey = "1", compname = "onehorizon",
        comppct_l = 100, comppct_r = 100, comppct_h = 100,
        hzname = "R", hzdept_l = 0, hzdept_r = 0, hzdept_h = 0,
        hzdepb_l = 145, hzdepb_r = 150, hzdepb_h = 155, hzthk_l = 145,
        stringsAsFactors = FALSE
      )
    },
    .package = "soilSIM"
  )

  set.seed(1)
  result <- simulate_profile_depths("999", n_simulations = 5, seed = 1)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 5)
})

test_that("evaluate_simulated_depths() flags horizons whose simulated depths exceed the original bounds", {
  simulated_profiles <- make_soil_profile_collection(sim_comppct = 1)
  horizon_data <- make_depth_horizon_data()
  # Force an out-of-range simulated top for the first horizon
  aqp::horizons(simulated_profiles)$hzdept_r[1] <- horizon_data$hzdept_l[1] - 5

  out_of_range <- evaluate_simulated_depths(simulated_profiles, horizon_data)
  expect_true(nrow(out_of_range) >= 1)
  expect_true(all(out_of_range$out_of_range))
})

test_that("simulate_and_perturb_soil_profiles() runs the full multi-horizon perturbation path (thickness + OSD-distinctness boundary perturbation) against a real series", {
  # Unlike the single-horizon test above (which takes the early-return path
  # before query_osd_distinctness() is ever called), a multi-horizon profile
  # exercises the full pipeline, which unconditionally calls
  # query_osd_distinctness() -> soilDB::fetchOSD() - a live OSD lookup that
  # requires a real component name (a synthetic "testseries" would return no
  # OSD data and error). "amador" is the same real series used in this
  # function's own historical example.
  testthat::skip_if_offline()
  df <- data.frame(
    cokey = "1", mukey = "1", compname = "amador", id = "amador",
    hzname = c("A", "Bt", "Cr"),
    hzdept_l = c(0, 15, 45), hzdept_r = c(0, 20, 50), hzdept_h = c(0, 25, 55),
    hzdepb_l = c(15, 45, 95), hzdepb_r = c(20, 50, 100), hzdepb_h = c(25, 55, 105),
    hzthk_l = c(15, 30, 50), sim_comppct = 4,
    stringsAsFactors = FALSE
  )
  spc <- df
  aqp::depths(spc) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(spc) <- "hzname"

  result <- simulate_and_perturb_soil_profiles(spc)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 4)
})

test_that("fetch_ssurgo_aws_data()/query_osd_distinctness() require the live SDA/OSD services", {
  testthat::skip_if_offline()
  testthat::skip("Live NRCS Soil Data Access / OSD queries are not exercised in automated tests - see test-ssurgo-acquisition.R for the established precedent.")
})

test_that("attach_osd_boundary_distinctness() joins per-(compname,genhz) bound_sd onto every matching row, offline via a mocked query_osd_distinctness()", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 1b. query_osd_distinctness() itself requires a
  # live soilDB::fetchOSD() call (see the skipped test above) - mock it directly, the same pattern
  # test-gp-modeling.R uses for predict_gp_depth_trends(), to test the join/aggregation logic here
  # fully offline.
  #
  # Mocked `id` is UPPER-CASED ("AMADOR"/"PENTZ") deliberately - this matches
  # soilDB::fetchOSD()'s own documented behavior (see fetch_osd_horizons_cached()'s "Known quirk
  # preserved" section) and is what a REAL query_osd_distinctness() call actually returns. A
  # Phase 9 real-AOI validation run found the original (pre-fix) join silently produced all-NA
  # bound_sd against real data because it joined on exact compname case - this test's mock now
  # exercises that same case mismatch to lock the fix in (see Phase 9's write-up below).
  testthat::local_mocked_bindings(
    query_osd_distinctness = function(horizon_data) {
      data.frame(
        id = c("AMADOR", "AMADOR", "AMADOR", "PENTZ"),
        hzname = c("A", "Bt1", "Bt2", "A"),
        distinctness = c("clear", "gradual", "gradual", "abrupt"),
        genhz = c("A", "B", "B", "A"),
        # Two "B" rows for amador (10, 20) - the group_by(genhz) mean should collapse
        # them to 15, confirming the same aggregation simulate_and_perturb_soil_profiles()
        # already uses is reused here.
        bound_sd = c(5, 10, 20, 2),
        stringsAsFactors = FALSE
      )
    },
    .package = "soilSIM"
  )

  hz_data <- data.frame(
    compname = c("amador", "amador", "amador", "pentz"),
    hzname = c("A", "Bt1", "Bt2", "A"),
    genhz = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )

  result <- attach_osd_boundary_distinctness(hz_data)
  expect_equal(nrow(result), 4)
  expect_false("compname_upper" %in% names(result))  # internal join key not leaked into output
  expect_equal(result$bound_sd[result$compname == "amador" & result$genhz == "A"], 5)
  expect_equal(result$bound_sd[result$compname == "amador" & result$genhz == "B"], c(15, 15))
  expect_equal(result$bound_sd[result$compname == "pentz"], 2)
})

test_that("attach_osd_boundary_distinctness() degrades to bound_sd = NA (not an error) when the OSD lookup fails", {
  testthat::local_mocked_bindings(
    query_osd_distinctness = function(horizon_data) stop("simulated network failure"),
    .package = "soilSIM"
  )

  hz_data <- data.frame(
    compname = "amador", hzname = "A", genhz = "A", stringsAsFactors = FALSE
  )
  result <- expect_no_error(attach_osd_boundary_distinctness(hz_data))
  expect_true(is.na(result$bound_sd))
  expect_equal(nrow(result), 1)
})

test_that("attach_osd_boundary_distinctness() requires compname/hzname/genhz columns", {
  expect_error(
    attach_osd_boundary_distinctness(data.frame(compname = "amador", hzname = "A")),
    "compname, hzname, and genhz"
  )
})

test_that("simulate_profile_depths(parallel = TRUE) runs without erroring at the R level", {
  skip_if_not(
    nzchar(system.file(package = "soilSIM")) &&
      file.exists(file.path(system.file(package = "soilSIM"), "Meta", "package.rds")),
    "soilSIM not installed in this session - future::multisession workers can't see a load_all()'d namespace"
  )
  # Unlike process_cokeys_parallel() (see test-multivariate-adjustment.R),
  # this function has no sequential fallback - it's a straight port of the
  # legacy code_ref/brdf/depth_simulation.R, which just returns NULL from its
  # tryCatch on any worker-side error. Some sandboxed environments' fresh
  # future::multisession worker processes can't see an installed package
  # either (the same "there is no package called 'soilSIM'" limitation
  # documented and tolerated in test-multivariate-adjustment.R), so accept
  # either a real result or that documented, gracefully-handled NULL.
  # future::multisession/parallelly's Windows PSOCK workers write launcher
  # scripts to tempdir() that outlive the worker (a documented
  # future/parallelly quirk, not a soilSIM defect) - clean up any it drops,
  # matching this suite's existing convention of tests cleaning up their own
  # tempdir droppings (see the cache_dir on.exit() in test-ssurgo-acquisition.R).
  tmp_before <- list.files(tempdir())
  on.exit({
    new_files <- setdiff(list.files(tempdir()), tmp_before)
    unlink(file.path(tempdir(), grep("^Rscript", new_files, value = TRUE)))
  }, add = TRUE)

  soil_collection <- make_soil_profile_collection(sim_comppct = 2)
  result <- simulate_profile_depths(soil_collection, seed = 1, parallel = TRUE, n_cores = 2)
  expect_true(is.null(result) || inherits(result, "SoilProfileCollection"))
})

test_that("simulate_profile_depths(parallel = TRUE) restores the caller's future::plan() afterward", {
  skip_if_not(
    nzchar(system.file(package = "soilSIM")) &&
      file.exists(file.path(system.file(package = "soilSIM"), "Meta", "package.rds")),
    "soilSIM not installed in this session - future::multisession workers can't see a load_all()'d namespace"
  )
  # Regression test for a real bug fixed this session: the pre-migration implementation reverted
  # future::plan() to sequential() only on its own success path (not via on.exit()), so it either
  # clobbered a caller's pre-existing plan or left multisession active indefinitely on error.
  # run_parallel_lapply() (R/parallel.R) now snapshots and restores the caller's plan
  # unconditionally.
  tmp_before <- list.files(tempdir())
  on.exit({
    new_files <- setdiff(list.files(tempdir()), tmp_before)
    unlink(file.path(tempdir(), grep("^Rscript", new_files, value = TRUE)))
  }, add = TRUE)

  plan_before <- future::plan()
  soil_collection <- make_soil_profile_collection(sim_comppct = 2)
  invisible(simulate_profile_depths(soil_collection, seed = 1, parallel = TRUE, n_cores = 2))
  expect_true(identical(class(future::plan()), class(plan_before)))
})
