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

test_that("sim_component_comp() adds a sim_comppct column, one row per distinct component", {
  data <- make_component_data()
  result <- sim_component_comp(data, n_simulations = 500)
  expect_equal(nrow(result), 2)
  expect_true("sim_comppct" %in% names(result))
  expect_true(all(result$sim_comppct > 0))
  # sim_comppct ~= round(n_simulations * comppct_r / 100) in expectation
  expect_equal(result$sim_comppct[1], round(500 * 55 / 100), tolerance = 20)
})

test_that("sim_component_comp() fills missing comppct_l/comppct_h from comppct_r +/- 2", {
  data <- data.frame(
    mukey = "1", cokey = "1", compname = "compA",
    comppct_l = NA_real_, comppct_r = 50, comppct_h = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- sim_component_comp(data, n_simulations = 200)
  expect_equal(nrow(result), 1)
  expect_true(result$sim_comppct > 0)
})

test_that("sim_component_comp() output joins onto horizon data by cokey (documented grain-mismatch fix)", {
  component_data <- make_component_data()
  comp_result <- sim_component_comp(component_data, n_simulations = 100)

  horizon_data <- data.frame(
    cokey = c("1", "1", "2"), hzname = c("A", "Bt", "A"),
    hzdept_r = c(0, 20, 0), hzdepb_r = c(20, 50, 30),
    stringsAsFactors = FALSE
  )
  joined <- dplyr::left_join(horizon_data, comp_result[, c("cokey", "sim_comppct")], by = "cokey")
  expect_true(all(!is.na(joined$sim_comppct)))
  expect_equal(nrow(joined), nrow(horizon_data))
})

test_that("calculate_mode() returns the most frequent value", {
  expect_equal(calculate_mode(c(1, 2, 2, 3, 3, 3, 4)), 3)
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
          list(ilr1_lrh = c(min(sim_txt_ilr[, "z1"]), calculate_mode(sim_txt_ilr[, "z1"]), max(sim_txt_ilr[, "z1"])),
               ilr2_lrh = c(min(sim_txt_ilr[, "z2"]), calculate_mode(sim_txt_ilr[, "z2"]), max(sim_txt_ilr[, "z2"])))
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
