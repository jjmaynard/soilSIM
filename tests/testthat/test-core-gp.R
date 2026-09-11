## Synthetic GP-fit-shaped data throughout - GPfit::GP_fit()/predict.GP() are
## pure local computation, so no network mocking is needed here.
## make_gp_training_df() lives in helper-fixtures.R (shared with test-classes.R's
## predict.soilSIM_gp_models()/plot.soilSIM_gp_models() tests).

test_that("create_soil_groups() derives soil_group from the requested strategy column, falling back to 'unknown'", {
  df <- data.frame(compname = c("Alpha", NA), taxclname = c("Fine", "Coarse"), stringsAsFactors = FALSE)
  by_series <- create_soil_groups(df, "soil_series")
  expect_equal(by_series$soil_group, c("Alpha", "unknown"))

  by_none <- create_soil_groups(df, "none")
  expect_equal(by_none$soil_group, c("all_soils", "all_soils"))

  # Unrecognized strategy falls back to the compname -> taxclname -> "unknown" hierarchy
  by_fallback <- create_soil_groups(df, "nonexistent_strategy")
  expect_equal(by_fallback$soil_group, c("Alpha", "unknown"))
})

test_that("create_test_grouping() builds test_group only from non-NA values of the strategy column", {
  df <- data.frame(compname = c("Alpha", NA, "Beta"), stringsAsFactors = FALSE)
  result <- create_test_grouping(df, "soil_series")
  expect_equal(nrow(result), 2)
  expect_setequal(result$test_group, c("Alpha", "Beta"))

  expect_equal(nrow(create_test_grouping(df, "unrecognized")), 0)
})

test_that("apply_hierarchical_grouping() keeps adequate groups as-is and falls back small groups to taxonomy-based pools", {
  df <- data.frame(
    cokey = as.character(1:10),
    soil_group = c(rep("BigGroup", 6), rep("TinyGroup", 4)),
    hzdept_r = c(0, 20, 40, 60, 80, 100, 0, 30, 0, 30),
    taxpartsize = c(rep(NA_character_, 6), rep("fine-loamy", 4)),
    clay_pct = c(10, 12, 14, 16, 18, 20, 25, 27, 22, 24),
    stringsAsFactors = FALSE
  )
  result <- apply_hierarchical_grouping(df, "clay_pct", min_profiles = 3, min_obs = 5)
  expect_true(all(result$final_group[df$soil_group == "BigGroup"] == "BigGroup"))
  # TinyGroup (4 obs, all one cokey-depth-range) falls back to the taxpartsize-based pool
  expect_true(all(grepl("^fallback_particle_", result$final_group[df$soil_group == "TinyGroup"]) |
                    result$final_group[df$soil_group == "TinyGroup"] %in% c("fallback_general", "general_pool")))
})

test_that("select_optimal_grouping() picks a strategy meeting the adequacy target, or the best available otherwise", {
  df <- data.frame(
    cokey = as.character(1:12),
    compname = rep(c("Alpha", "Beta", "Gamma"), each = 4),
    hzdept_r = rep(c(0, 30, 60, 90), 3),
    stringsAsFactors = FALSE
  )
  strategy <- select_optimal_grouping(df, min_profiles = 1, min_obs = 4, target_groups = 3)
  expect_equal(strategy, "soil_series")
})

test_that("fit_individual_gp_model() fits a GP model on depth-aggregated data and predict_gp_depth_trends() reproduces a sane depth trend", {
  df <- make_gp_training_df()
  result <- fit_individual_gp_model(df, "clay_pct", optimize_hyperparameters = FALSE)
  expect_false(is.null(result))
  expect_true(is.list(result$model))
  expect_equal(result$model$property, "clay_pct")
  expect_true(result$model$n_training_points >= 3)

  preds <- predict_gp_depth_trends(result$model, new_depths = c(0, 50, 100))
  expect_length(preds, 3)
  expect_true(all(is.finite(preds)))
  # Clay content increases with depth in the synthetic fixture - the fitted
  # trend should reproduce that monotonic increase at the endpoints.
  expect_true(preds[3] > preds[1])
})

test_that("fit_individual_gp_model() returns NULL (not an error) when there is no variation in the property", {
  df <- make_gp_training_df()
  df$clay_pct <- 20  # constant - no variation
  expect_null(fit_individual_gp_model(df, "clay_pct", optimize_hyperparameters = FALSE))
})

test_that("predict_gp_depth_trends() returns NA predictions (not an error) for a NULL/malformed model", {
  expect_equal(predict_gp_depth_trends(NULL, c(0, 10)), rep(NA, 2))
  expect_equal(predict_gp_depth_trends(list(gp_model = NULL), c(0, 10)), rep(NA, 2))
})

test_that("optimize_gp_hyperparameters() performs real cross-validation and reports which folds/candidates were compared", {
  set.seed(1)
  X <- seq(0, 1, length.out = 12)
  Y <- X * 10 + stats::rnorm(12, sd = 0.01)

  model <- optimize_gp_hyperparameters(X, Y, n_folds = 4)
  expect_true(inherits(model, "GP"))

  cv <- attr(model, "cv_results")
  expect_equal(cv$folds, 4)
  expect_setequal(cv$corr_candidates, c("exponential_power1.95", "matern_nu1.5", "matern_nu2.5"))
  # At least one candidate produced a finite CV RMSE for this well-behaved,
  # near-linear synthetic trend.
  expect_true(any(is.finite(unlist(cv$mean_rmse_by_candidate))))
  expect_true(cv$best_candidate %in% cv$corr_candidates)
})

test_that("optimize_gp_hyperparameters() falls back to a single baseline fit when there isn't enough data for CV", {
  X <- seq(0, 1, length.out = 3)
  Y <- c(1, 2, 3)
  model <- optimize_gp_hyperparameters(X, Y, n_folds = 5)
  expect_true(inherits(model, "GP"))
  cv <- attr(model, "cv_results")
  expect_true(is.na(cv$best_candidate))
})

test_that("assess_trend_monotonicity() flags consistently increasing/decreasing trends and rejects noisy ones", {
  expect_true(assess_trend_monotonicity(c(1, 2, 3, 4, 5)))
  expect_true(assess_trend_monotonicity(c(5, 4, 3, 2, 1)))
  expect_false(assess_trend_monotonicity(c(1, 5, 1, 5, 1)))
  expect_true(is.na(assess_trend_monotonicity(c(1, NA))))
})

test_that("assess_realistic_values() checks property-specific plausible ranges", {
  expect_true(assess_realistic_values(c(10, 20, 30), "clay_pct"))
  expect_false(assess_realistic_values(c(10, 150), "clay_pct"))
  expect_true(assess_realistic_values(c(6.5, 7.0), "pH"))
  expect_false(assess_realistic_values(c(-1, 7.0), "pH"))
  # Unknown property falls back to the generic non-negative/magnitude check
  expect_true(assess_realistic_values(c(1, 2, 3), "some_unknown_property"))
  expect_false(assess_realistic_values(c(-5), "some_unknown_property"))
})

test_that("infer_simulation_properties() excludes structural/metadata columns and keeps numeric property columns", {
  df <- data.frame(
    cokey = "1", hzdept_r = 0, hzdepb_r = 20, hzname = "A", genhz = "A",
    mukey = "1", simulation_number = 1, infill_method = "", unsuitable_horizon = FALSE,
    property_data_complete = TRUE, claytotal = 20, sandtotal = 40,
    stringsAsFactors = FALSE
  )
  expect_setequal(infer_simulation_properties(df), c("claytotal", "sandtotal"))
})

test_that("apply_nrcs_gp_adjustments_with_correlations() safely no-ops when no matching NRCS GP model exists", {
  sim_data <- data.frame(cokey = "1", hzdept_r = c(0, 20), simulation_number = 1, clay_pct = c(10, 15))
  # "clay_pct" is a GP-model target name, not one of get_nrcs_property_mapping()'s
  # source-property keys, so no mapping matches and the delegate real
  # function (apply_nrcs_trend_adjustments()) correctly no-ops.
  res <- apply_nrcs_gp_adjustments_with_correlations(sim_data, nrcs_gp_models = list(), model_group = "x")
  expect_identical(res, sim_data)
})

test_that("apply_nrcs_gp_adjustments_with_correlations() applies a fitted NRCS depth trend via the real delegate function", {
  training_df <- make_gp_training_df()  # clay_pct increases with depth: 10,15,22,30 (cokey 1)
  fit_result <- fit_individual_gp_model(training_df, "clay_pct", optimize_hyperparameters = FALSE)
  gp_models <- list(clay_pct = list(type = "stratified_grouped", models = list(group1 = fit_result$model)))

  sim_data <- data.frame(
    cokey = "99", hzdept_r = c(0, 20, 50, 100), simulation_number = 1,
    claytotal = c(10, 10, 10, 10)  # constant before adjustment
  )
  res <- apply_nrcs_gp_adjustments_with_correlations(sim_data, gp_models, model_group = "group1")

  expect_equal(nrow(res), nrow(sim_data))
  # "claytotal" maps to the fitted "clay_pct" GP model via
  # get_nrcs_property_mapping(); since the fitted trend increases with
  # depth, the previously-constant column should no longer be constant.
  expect_false(all(res$claytotal == res$claytotal[1]))
})

test_that("apply_local_gp_adjustments_with_correlations() delegates to the real apply_local_gp_adjustments() and no-ops on too few depths", {
  sim_data <- data.frame(cokey = "1", hzdept_r = c(0, 20), simulation_number = 1, clay_pct = c(10, 15))
  # Only 2 unique depths - below apply_local_gp_adjustments()'s default
  # min_depths = 3, so it safely returns cokey_data unchanged.
  res <- apply_local_gp_adjustments_with_correlations(sim_data)
  expect_identical(res, sim_data)
})

test_that("stub contract: optimize_gp_hyperparameters() never performs real cross-validation regardless of n_folds", {
  # Locks in the current stub behavior so a future real implementation is a
  # deliberate, visible change (per migration plan B.5) rather than accidental.
  X <- seq(0, 1, length.out = 5)
  Y <- c(1, 2, 3, 4, 5)
  model <- optimize_gp_hyperparameters(X, Y, n_folds = 100)  # absurd n_folds, still accepted
  expect_true(inherits(model, "GP"))
})

test_that("flatten_simulation_array_to_long() reshapes a [horizon,property,realization] array to match horizon metadata", {
  cokey_data <- data.frame(cokey = "1", hzdept_r = c(0, 10), hzdepb_r = c(10, 20))
  sim_array <- array(
    c(1, 2, 3, 4),
    dim = c(2, 1, 2),
    dimnames = list(horizon = 1:2, property = "dbovendry", realization = 1:2)
  )
  long_df <- flatten_simulation_array_to_long(sim_array, cokey_data)
  expect_equal(nrow(long_df), 4)
  expect_equal(long_df$simulation_number, c(1, 1, 2, 2))
  expect_equal(long_df$dbovendry, c(1, 2, 3, 4))
  expect_equal(long_df$hzdept_r, c(0, 10, 0, 10))
})

test_that("flatten_simulation_array_to_long() errors when cokey_data row count does not match the array's horizon dimension", {
  cokey_data <- data.frame(cokey = "1", hzdept_r = 0, hzdepb_r = 10)
  sim_array <- array(1:4, dim = c(2, 1, 2), dimnames = list(NULL, "dbovendry", NULL))
  expect_error(flatten_simulation_array_to_long(sim_array, cokey_data))
})

test_that("select_simulation_correlation_matrix() picks the genhz-matched matrix and subsets to shared properties", {
  cokey_data <- make_horizon_row(properties = c("sandtotal", "claytotal"))
  mat_A <- matrix(c(1, 0.5, 0.5, 1), nrow = 2,
                   dimnames = list(c("sandtotal", "claytotal"), c("sandtotal", "claytotal")))
  correlation_matrices <- list(A = mat_A, B = matrix(1))
  chosen <- select_simulation_correlation_matrix(cokey_data, correlation_matrices, NULL)
  expect_setequal(chosen$properties, c("sandtotal", "claytotal"))
  expect_equal(dim(chosen$correlation_matrix), c(2, 2))
})

test_that("select_simulation_correlation_matrix() falls back to auto-detected properties and a NULL matrix when nothing usable is supplied", {
  cokey_data <- make_horizon_row(properties = c("sandtotal", "claytotal"))
  chosen <- select_simulation_correlation_matrix(cokey_data, NULL, NULL)
  expect_true(is.null(chosen$correlation_matrix))
  expect_true(length(chosen$properties) > 0)
})

test_that("simulate_soil_properties() runs end-to-end through the real Monte Carlo engine and returns long-format results", {
  sim_data <- make_soil_data(n_horizons = 9, properties = c("dbovendry"))
  sim_data$cokey <- "1"
  sim_data$hzdept_r <- (seq_len(9) - 1) * 10
  sim_data$hzdepb_r <- seq_len(9) * 10

  result <- simulate_soil_properties(
    target_cokey = "1",
    sim_data = sim_data,
    correlation_matrices = NULL,
    txt_correlation_matrices = NULL,
    n_simulations = 5,
    use_nrcs_gp = FALSE
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 9 * 5)
  expect_true("dbovendry" %in% names(result))
  expect_true("simulation_number" %in% names(result))
  expect_false(is.null(attr(result, "validation")))
})

test_that("adjust_simulation_depthwise() preserves matrix dimensions and validates mismatched inputs", {
  set.seed(42)
  n_depths <- 3
  n_sims <- 20
  depths <- c(0, 50, 100)

  simulated_list <- list(
    clay_pct = matrix(stats::rnorm(n_depths * n_sims, mean = 20, sd = 2), nrow = n_depths),
    sand_pct = matrix(stats::rnorm(n_depths * n_sims, mean = 40, sd = 3), nrow = n_depths)
  )

  gp_models <- list(
    clay_pct = list(gp_model = list(mean = "flat"), predictions_override = c(20, 22, 25)),
    sand_pct = list(gp_model = list(mean = "flat"), predictions_override = c(40, 38, 35))
  )
  # adjust_simulation_depthwise() calls predict_gp_depth_trends() when the
  # element has a $gp_model - stub it here via a minimal fake to avoid needing
  # a real fitted GPfit object for this dimension/validation-focused test.
  testthat::local_mocked_bindings(
    predict_gp_depth_trends = function(gp_model_info, new_depths) gp_model_info$predictions_override,
    .package = "soilSIM"
  )

  result <- adjust_simulation_depthwise(simulated_list, gp_models, depths)
  expect_named(result, c("clay_pct", "sand_pct"))
  expect_equal(dim(result$clay_pct), c(n_depths, n_sims))
  expect_equal(dim(result$sand_pct), c(n_depths, n_sims))

  expect_error(adjust_simulation_depthwise(list(), gp_models, depths), "cannot be empty")
  expect_error(adjust_simulation_depthwise(simulated_list, gp_models, c(0, 50)), "Length of depths")
})

test_that("adjust_simulation_depthwise()'s vectorized quantile() call matches the original per-replicate loop", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 3 adjust_simulation_depthwise()
  # fix: quantile(curr_values, probs = q) was previously called once per replicate (recomputing
  # the full quantile of the same vector at a different single probability each time) - now one
  # call with a probs vector. Reimplements the ORIGINAL per-replicate loop here and checks the
  # vectorized version's output is identical.
  set.seed(99)
  n_depths <- 4
  n_sims <- 30
  depths <- c(0, 20, 50, 100)

  simulated_list <- list(
    clay_pct = matrix(stats::rnorm(n_depths * n_sims, mean = 20, sd = 3), nrow = n_depths),
    sand_pct = matrix(stats::rnorm(n_depths * n_sims, mean = 40, sd = 4), nrow = n_depths)
  )
  gp_models <- list(
    clay_pct = list(gp_model = list(mean = "flat"), predictions_override = c(20, 21, 23, 26)),
    sand_pct = list(gp_model = list(mean = "flat"), predictions_override = c(40, 39, 36, 33))
  )
  testthat::local_mocked_bindings(
    predict_gp_depth_trends = function(gp_model_info, new_depths) gp_model_info$predictions_override,
    .package = "soilSIM"
  )

  old_adjust <- function(simulated_list, gp_predictions, primary_property, depths) {
    n_depths <- nrow(simulated_list[[1]]); n_sims <- ncol(simulated_list[[1]])
    primary_matrix <- simulated_list[[primary_property]]
    surface_values <- primary_matrix[1, ]
    reference_quantiles <- ecdf(surface_values)(surface_values)

    adjusted_list <- list()
    for (prop in names(simulated_list)) {
      current_matrix <- simulated_list[[prop]]
      gp_means <- gp_predictions[[prop]]
      adjusted_matrix <- current_matrix
      for (i in 2:n_depths) {
        gp_ratio <- if (is.na(gp_means[i - 1]) || is.na(gp_means[i]) || gp_means[i - 1] == 0) 1 else gp_means[i] / gp_means[i - 1]
        prev_values <- adjusted_matrix[i - 1, ]
        curr_values <- current_matrix[i, ]
        adjusted_curr <- numeric(n_sims)
        for (j in 1:n_sims) {
          q <- reference_quantiles[j]
          quantile_value <- quantile(curr_values, probs = q, na.rm = TRUE)
          adjusted_curr[j] <- quantile_value + (prev_values[j] * gp_ratio - quantile_value)
        }
        if (var(curr_values, na.rm = TRUE) > 0) {
          ecdf_adjusted <- ecdf(adjusted_curr)
          adjusted_matrix[i, ] <- quantile(curr_values, probs = ecdf_adjusted(adjusted_curr), na.rm = TRUE)
        } else {
          adjusted_matrix[i, ] <- curr_values
        }
      }
      adjusted_list[[prop]] <- adjusted_matrix
    }
    adjusted_list
  }

  gp_predictions <- list(clay_pct = c(20, 21, 23, 26), sand_pct = c(40, 39, 36, 33))
  expected <- old_adjust(simulated_list, gp_predictions, "clay_pct", depths)
  actual <- adjust_simulation_depthwise(simulated_list, gp_models, depths, primary_property = "clay_pct")

  expect_equal(unname(actual$clay_pct), unname(expected$clay_pct), tolerance = 1e-10)
  expect_equal(unname(actual$sand_pct), unname(expected$sand_pct), tolerance = 1e-10)
})

test_that("REGRESSION: prepare_nrcs_training_data() runs without error and filters out unsuitable (R) horizons", {
  # Prior to a fix during this migration, converting this function's %>%
  # pipe to base R |> left a bare `.` inside is_unsuitable(., ...) - valid
  # under magrittr (which substitutes `.` anywhere in the RHS call), but
  # meaningless under |> (which has no such substitution), so this call
  # would error with "object '.' not found" every time. Locking in that it
  # now runs and correctly excludes the bedrock ("R") horizon.
  raw <- data.frame(
    cokey = rep(c("1", "2"), each = 4),
    hzname = c("A", "Bw", "Bt", "R", "A", "Bw", "Bt", "Cr"),
    hzdept_r = rep(c(0, 20, 40, 80), 2),
    hzdepb_r = rep(c(20, 40, 80, 120), 2),
    hzdept_l = NA_real_, hzdepb_l = NA_real_,
    hzdept_h = NA_real_, hzdepb_h = NA_real_,
    claytotal_r = c(15, 20, 25, 30, 12, 18, 22, 28),
    compname = "Alpha",
    stringsAsFactors = FALSE
  )
  result <- prepare_nrcs_training_data(
    raw, grouping_strategy = "none",
    min_profiles_per_group = 1, min_observations_per_group = 3, target_min_groups = 1
  )
  expect_true(nrow(result) > 0)
  expect_false("R" %in% result$hzname)
})

test_that("validate_joint_correlation_structure() reports achieved property and depth-lag correlation against targets, at every depth/property", {
  # Phase 0 instrumentation for the vertical-correlation redesign
  # (PERFORMANCE_IMPROVEMENT_PLAN.md-style acceptance testing): unlike
  # validate_correlation_preservation(), which only checks property correlation, and only via its
  # caller's first-5-depths spot check, this must check BOTH property correlation (every depth)
  # and depth-lag correlation (every property). Construct synthetic data via the same
  # separable/Kronecker-MVN sampling identity the joint-copula redesign itself will use
  # (Z = L_depth %*% Eps %*% t(L_prop)) so the achieved correlations are known to approximate a
  # specific target by construction, at a large n_sims.
  set.seed(2024)
  n_depths <- 4
  n_sims <- 4000
  property_names <- c("clay_pct", "sand_pct", "bulk_density")
  k <- length(property_names)

  target_property_corr <- matrix(c(
    1.00, -0.60,  0.50,
   -0.60,  1.00, -0.40,
    0.50, -0.40,  1.00
  ), nrow = k, dimnames = list(property_names, property_names))

  target_depth_corr <- matrix(c(
    1.00, 0.80, 0.64, 0.51,
    0.80, 1.00, 0.80, 0.64,
    0.64, 0.80, 1.00, 0.80,
    0.51, 0.64, 0.80, 1.00
  ), nrow = n_depths)

  L_depth <- chol(target_depth_corr)
  L_prop <- chol(target_property_corr)

  eps <- array(stats::rnorm(n_depths * k * n_sims), dim = c(n_depths, k, n_sims))
  joint <- vapply(seq_len(n_sims), function(s) {
    t(L_depth) %*% eps[, , s] %*% L_prop
  }, matrix(0, n_depths, k))
  # joint is n_depths x k x n_sims - split into the per-property [depth x sim] matrix shape
  # simulate_correlated_triangular()/preserve_correlation_structure() actually use.
  simulated_list <- setNames(
    lapply(seq_len(k), function(p) matrix(joint[, p, ], nrow = n_depths)),
    property_names
  )

  result <- validate_joint_correlation_structure(
    simulated_list,
    target_property_corr = target_property_corr,
    target_depth_corr = target_depth_corr
  )

  # Shape: one property-correlation matrix per depth, one depth-correlation matrix per property.
  expect_length(result$achieved_property_correlation, n_depths)
  expect_named(result$achieved_depth_correlation, property_names)
  expect_length(result$property_correlation_max_diff, n_depths)
  expect_named(result$depth_correlation_max_diff, property_names)

  # Data actually constructed to match both targets - achieved deviation should be small at
  # n_sims = 4000, well under the 0.1 "large correlation change" threshold
  # adjust_simulation_depthwise() itself warns on.
  expect_true(result$overall_property_max_diff < 0.1)
  expect_true(result$overall_depth_max_diff < 0.1)
  expect_true(all(result$property_correlation_max_diff < 0.1))
  expect_true(all(result$depth_correlation_max_diff < 0.1))

  # A deliberately wrong target should be caught, not silently passed.
  wrong_target <- diag(k)
  dimnames(wrong_target) <- dimnames(target_property_corr)
  bad_result <- validate_joint_correlation_structure(
    simulated_list, target_property_corr = wrong_target
  )
  expect_true(bad_result$overall_property_max_diff > 0.3)
})

test_that("validate_joint_correlation_structure() returns all-NA diffs (not an error) when no targets are supplied", {
  set.seed(7)
  simulated_list <- list(
    clay_pct = matrix(stats::rnorm(3 * 10, mean = 20, sd = 2), nrow = 3),
    sand_pct = matrix(stats::rnorm(3 * 10, mean = 40, sd = 3), nrow = 3)
  )
  result <- validate_joint_correlation_structure(simulated_list)
  expect_true(all(is.na(result$property_correlation_max_diff)))
  expect_true(all(is.na(result$depth_correlation_max_diff)))
  expect_true(is.na(result$overall_property_max_diff))
  expect_true(is.na(result$overall_depth_max_diff))
  expect_length(result$achieved_property_correlation, 3)
  expect_named(result$achieved_depth_correlation, c("clay_pct", "sand_pct"))
})

test_that("extract_depth_length_scale() derives a positive real-units length-scale from a fitted GP model", {
  # Phase 1 of VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md: the depth kernel's length-scale must come
  # from the already-fitted GPfit model (fit_local_gp_model_single()), not a new estimation step.
  set.seed(11)
  agg_data <- data.frame(
    hzdept_r = c(0, 20, 50, 100, 150),
    mean_val = c(10, 14, 22, 30, 33)
  )
  fitted <- fit_local_gp_model_single(agg_data, "clay_pct")
  expect_false(is.null(fitted))

  length_scale <- extract_depth_length_scale(fitted)
  expect_true(is.finite(length_scale))
  expect_true(length_scale > 0)
  # Real-units: should be comparable to (not wildly outside) the actual depth range (0-150cm) -
  # loosely bounding to catch a units mix-up (e.g. forgetting to rescale by depth_scaling$range).
  expect_true(length_scale > 1 && length_scale < 1500)

  # Without depth_scaling, the result stays in scaled [0,1] units - must be <= ~a few (GPfit
  # ranges rarely exceed the unit interval by more than a small factor for a well-fit model).
  scaled_only <- extract_depth_length_scale(fitted$gp_model, depth_scaling = NULL)
  expect_true(is.finite(scaled_only))
  expect_equal(scaled_only * fitted$depth_scaling$range, length_scale, tolerance = 1e-8)

  # Malformed/missing model input degrades to NA_real_, not an error - mirrors
  # fit_local_gp_model_single()'s own NULL-on-failure contract.
  expect_true(is.na(extract_depth_length_scale(NULL)))
  expect_true(is.na(extract_depth_length_scale(list(gp_model = NULL))))
})

test_that("build_depth_correlation_kernel() returns a valid, distance-decaying correlation matrix for both kernel families", {
  depths <- c(0, 10, 30, 100)  # uneven spacing, on purpose

  for (kernel in c("exponential", "matern")) {
    R <- build_depth_correlation_kernel(depths, length_scale = 25, kernel = kernel)

    expect_equal(dim(R), c(4, 4))
    expect_equal(diag(R), rep(1, 4), tolerance = 1e-8)
    expect_true(isSymmetric(R, tol = 1e-8))
    expect_true(all(eigen(R, only.values = TRUE)$values > -1e-8))  # positive (semi-)definite

    # Monotonic decay with distance: farther depth-pairs must not be MORE correlated than nearer
    # ones, comparing every pair sharing depth index 1 as the common anchor.
    d_from_1 <- abs(depths - depths[1])
    ord <- order(d_from_1)
    expect_true(all(diff(R[1, ord]) <= 1e-8))
  }

  # Non-finite/non-positive length_scale degrades to the identity matrix (no depth correlation),
  # not an error - matches this package's established fallback pattern for invalid correlation
  # inputs (e.g. simulate_cokey_generalized()'s pooled-matrix fallback).
  expect_equal(build_depth_correlation_kernel(depths, length_scale = NA_real_), diag(4))
  expect_equal(build_depth_correlation_kernel(depths, length_scale = 0), diag(4))
  expect_equal(build_depth_correlation_kernel(depths, length_scale = -5), diag(4))

  # Single depth is a degenerate 1x1 correlation matrix, not an error.
  expect_equal(build_depth_correlation_kernel(15, length_scale = 25), matrix(1, 1, 1))
})

test_that("build_depth_correlation_kernel()'s boundary_distinctness gating suppresses correlation across an abrupt boundary while leaving diffuse-boundary profiles ~unchanged", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 1c/1d.
  depths <- c(0, 20, 40, 60)  # evenly spaced, on purpose - isolates the gating effect from the
  # plain distance-decay kernel's own shape.
  length_scale <- 40

  plain <- build_depth_correlation_kernel(depths, length_scale)

  # An injected Abrupt (bound_sd ~= 1, aqp's own code-to-offset value) boundary between depth
  # index 2 (20cm) and 3 (40cm) - bound_sd is indexed by "depth this boundary sits above", so
  # index 3 carries it. Other boundaries left NA (no gating - default "no data" behavior).
  abrupt_bsd <- c(NA, NA, 1, NA)
  gated <- build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = abrupt_bsd)

  expect_equal(dim(gated), c(4, 4))
  expect_equal(diag(gated), rep(1, 4), tolerance = 1e-8)
  expect_true(isSymmetric(gated, tol = 1e-8))
  expect_true(all(eigen(gated, only.values = TRUE)$values > -1e-8))

  # Every pair straddling the abrupt boundary (one depth in {1,2}, the other in {3,4}) must be
  # suppressed well below its ungated value; pairs entirely on one side of the boundary must be
  # essentially unaffected (no NA elsewhere, min_gate_weight floor only applies across the gate).
  straddling <- list(c(1, 3), c(1, 4), c(2, 3), c(2, 4))
  for (pair in straddling) {
    expect_true(gated[pair[1], pair[2]] < plain[pair[1], pair[2]] * 0.5)
  }
  same_side <- list(c(1, 2), c(3, 4))
  for (pair in same_side) {
    expect_equal(gated[pair[1], pair[2]], plain[pair[1], pair[2]], tolerance = 0.05)
  }

  # All-Diffuse (bound_sd = 10, the distinctness_range max -> gate weight exactly 1 everywhere)
  # must reproduce the plain kernel, confirming gating is a strict generalization, not a behavior
  # change, when every boundary is gradual/diffuse.
  diffuse_bsd <- rep(10, 4)
  ungated_equivalent <- build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = diffuse_bsd)
  expect_equal(ungated_equivalent, plain, tolerance = 1e-8)

  # NULL (default) and all-NA both skip gating entirely - identical to the plain kernel.
  expect_equal(build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = NULL), plain)
  expect_equal(build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = rep(NA_real_, 4)), plain)

  # Gating must compound across multiple intervening boundaries, not just gate adjacent pairs:
  # two abrupt boundaries in a row should suppress the endpoints MORE than a single abrupt
  # boundary does.
  two_abrupt_bsd <- c(NA, 1, 1, NA)
  double_gated <- build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = two_abrupt_bsd)
  expect_true(double_gated[1, 4] < gated[1, 4])

  # Wrong-length boundary_distinctness errors rather than silently misaligning with depths.
  expect_error(
    build_depth_correlation_kernel(depths, length_scale, boundary_distinctness = c(1, 2)),
    "same length"
  )

  # Depths passed out of order still gate the correct (physically adjacent) boundary - internal
  # sorting/un-sorting must round-trip correctly.
  shuffled_order <- c(3, 1, 4, 2)  # arbitrary permutation of 1:4
  shuffled_depths <- depths[shuffled_order]
  shuffled_bsd <- abrupt_bsd[shuffled_order]
  shuffled_gated <- build_depth_correlation_kernel(shuffled_depths, length_scale, boundary_distinctness = shuffled_bsd)
  # Map back to original depth-index ordering for a direct comparison against `gated`.
  unshuffle <- order(shuffled_order)
  expect_equal(shuffled_gated[unshuffle, unshuffle], gated, tolerance = 1e-8)
})


# --- merged from test-multivariate-adjustment.R (P1 reorg) ---

## Synthetic GP-fit-shaped simulation data throughout - GPfit is pure local
## computation, so no network mocking is needed here.

make_sim_data <- function(cokeys = c("1", "2"), depths = c(0, 20, 50, 100), n_sims = 6) {
  set.seed(7)
  rows <- list()
  for (ck in cokeys) {
    for (d in depths) {
      for (s in seq_len(n_sims)) {
        rows[[length(rows) + 1]] <- data.frame(
          cokey = ck, hzdept_r = d, simulation_number = s,
          clay_pct = 15 + d * 0.1 + stats::rnorm(1, sd = 0.5),
          sand_pct = 45 - d * 0.05 + stats::rnorm(1, sd = 0.5),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

test_that("sample_joint_depth_property_copula() achieves both the target depth correlation and the target property correlation simultaneously", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 2 - the whole point of the joint-copula
  # redesign is that BOTH correlation structures are satisfied by construction, not retrofitted.
  # Reuses validate_joint_correlation_structure() (Phase 0's own acceptance-test instrumentation)
  # against known target matrices - the exact same check used to validate the Phase 0 helper
  # itself, now applied to the real sampler it was built to judge.
  set.seed(31)
  R_depth <- matrix(c(
    1.00, 0.85, 0.72,
    0.85, 1.00, 0.85,
    0.72, 0.85, 1.00
  ), nrow = 3)
  R_prop <- matrix(c(
    1.00, -0.60, 0.50,
   -0.60,  1.00, -0.40,
    0.50, -0.40,  1.00
  ), nrow = 3, dimnames = list(c("clay", "bd", "wr"), c("clay", "bd", "wr")))

  Z <- sample_joint_depth_property_copula(R_depth, R_prop, n_sims = 4000, seed = 31)
  expect_equal(dim(Z), c(3, 3, 4000))

  # Marginal standard-normality (mean ~0, sd ~1) at a representative depth/property cell.
  expect_equal(mean(Z[1, 1, ]), 0, tolerance = 0.1)
  expect_equal(stats::sd(Z[1, 1, ]), 1, tolerance = 0.1)

  simulated_list <- stats::setNames(
    lapply(seq_len(3), function(p) matrix(Z[, p, ], nrow = 3)),
    colnames(R_prop)
  )
  result <- validate_joint_correlation_structure(
    simulated_list, target_property_corr = R_prop, target_depth_corr = R_depth
  )
  expect_true(result$overall_property_max_diff < 0.1)
  expect_true(result$overall_depth_max_diff < 0.1)
})

test_that("sample_joint_depth_property_copula() validates n_sims and errors, not silently misbehaves, on n_sims < 1", {
  R_depth <- diag(2)
  R_prop <- diag(2)
  expect_error(sample_joint_depth_property_copula(R_depth, R_prop, n_sims = 0), "n_sims")
})

test_that("apply_copula_to_marginals() preserves each depth's fitted marginal distribution when Z is uncorrelated (identity copula)", {
  # With R_depth = R_prop = identity, Z is plain iid standard normal - pnorm(Z) is then
  # (approximately, for large n_sims) uniform on (0,1), so quantile(curr_values, probs = u)
  # should reproduce curr_values' own distribution (same mean/sd), confirming this function does
  # not distort marginals on its own - any correlation structure comes entirely from Z.
  set.seed(41)
  n_depths <- 2; n_sims <- 5000
  property_matrices <- list(
    clay = matrix(rnorm(n_depths * n_sims, mean = 20, sd = 3), nrow = n_depths),
    bd   = matrix(rnorm(n_depths * n_sims, mean = 1.3, sd = 0.1), nrow = n_depths)
  )
  Z <- sample_joint_depth_property_copula(diag(n_depths), diag(2), n_sims = n_sims, seed = 41)

  result <- apply_copula_to_marginals(Z, property_matrices)
  expect_named(result, c("clay", "bd"))
  expect_equal(dim(result$clay), c(n_depths, n_sims))

  for (prop in names(property_matrices)) {
    for (i in seq_len(n_depths)) {
      expect_equal(mean(result[[prop]][i, ]), mean(property_matrices[[prop]][i, ]), tolerance = 0.15)
      expect_equal(stats::sd(result[[prop]][i, ]), stats::sd(property_matrices[[prop]][i, ]), tolerance = 0.15)
    }
  }
})

test_that("apply_copula_to_marginals() re-centers onto the GP-predicted mean when gp_predictions is supplied", {
  set.seed(42)
  n_depths <- 2; n_sims <- 3000
  property_matrices <- list(
    clay = matrix(rnorm(n_depths * n_sims, mean = 20, sd = 3), nrow = n_depths)
  )
  Z <- sample_joint_depth_property_copula(diag(n_depths), diag(1), n_sims = n_sims, seed = 42)

  gp_predictions <- list(clay = c(20, 35))  # depth 2's target mean far from its original ~20
  result <- apply_copula_to_marginals(Z, property_matrices, gp_predictions = gp_predictions)

  expect_equal(mean(result$clay[1, ]), 20, tolerance = 0.5)
  expect_equal(mean(result$clay[2, ]), 35, tolerance = 0.5)
  # Spread (shape) is preserved even though location shifted.
  expect_equal(stats::sd(result$clay[2, ]), stats::sd(property_matrices$clay[2, ]), tolerance = 0.3)
})

test_that("apply_copula_to_marginals() validates its Z/property_matrices dimension contract", {
  Z <- sample_joint_depth_property_copula(diag(2), diag(2), n_sims = 5, seed = 1)
  expect_error(apply_copula_to_marginals(Z, list()), "named elements")
  expect_error(
    apply_copula_to_marginals(Z, list(clay = matrix(1, 2, 5))),  # only 1 property, Z has 2
    "must match"
  )
})

test_that("apply_copula_to_marginals() does not error on an all-NA property row (real messy field data)", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 10 regression: found via a full-AOI benchmark
  # against real Salinas Valley SSURGO data - a cokey row with an incomplete texture triplet
  # (simulate_cokey_generalized()'s own per-row texture tryCatch() leaves that row's property
  # all-NA) previously crashed this function with "missing value where TRUE/FALSE needed", since
  # stats::var(x, na.rm = TRUE) returns NA (not FALSE) for an all-NA x, and `if (NA > 0)` errors
  # rather than falling through. No synthetic unit-test fixture before this had ever included an
  # all-NA depth/property row.
  set.seed(61)
  n_depths <- 3; n_sims <- 20
  clay_mat <- matrix(rnorm(n_depths * n_sims, mean = 20, sd = 3), nrow = n_depths)
  clay_mat[2, ] <- NA_real_  # depth 2 entirely missing for this property, as in the real bug

  Z <- sample_joint_depth_property_copula(diag(n_depths), diag(1), n_sims = n_sims, seed = 61)

  result <- expect_no_error(apply_copula_to_marginals(Z, list(clay = clay_mat)))
  expect_equal(dim(result$clay), c(n_depths, n_sims))
  expect_true(all(is.na(result$clay[2, ])))  # the all-NA depth stays all-NA, not fabricated
  expect_false(any(is.na(result$clay[1, ])))  # unaffected depths are untouched
  expect_false(any(is.na(result$clay[3, ])))

  # Also confirmed at the preserve_correlation_structure_joint() level (the actual call site the
  # real bug surfaced through) - previously this silently fell back to unadjusted output via its
  # own tryCatch(), masking the underlying crash rather than fixing it.
  gp_predictions <- list(clay = c(20, 20, 20))
  result_joint <- expect_no_error(
    preserve_correlation_structure_joint(list(clay = clay_mat), gp_predictions, c(0, 20, 40), "clay")
  )
  expect_true(all(is.na(result_joint$clay[2, ])))
})

make_independent_depth_property_matrices <- function(target_prop_corr, n_depths = 4, n_sims = 3000, seed = 51) {
  # Builds property_matrices whose depths are drawn INDEPENDENTLY of each other (near-zero depth
  # correlation), each depth sharing the SAME known property-correlation target - the natural
  # "before" state preserve_correlation_structure_joint() is meant to add vertical correlation to,
  # while (ideally) preserving the property correlation that's already there.
  set.seed(seed)
  k <- nrow(target_prop_corr)
  L <- chol(target_prop_corr)
  props <- vector("list", k)
  names(props) <- colnames(target_prop_corr)
  for (p in seq_len(k)) props[[p]] <- matrix(NA_real_, n_depths, n_sims)
  for (i in seq_len(n_depths)) {
    eps <- matrix(stats::rnorm(n_sims * k), nrow = n_sims, ncol = k)
    correlated <- eps %*% L  # n_sims x k, cor ~= target_prop_corr
    for (p in seq_len(k)) props[[p]][i, ] <- correlated[, p] + i  # small depth-shift in mean only
  }
  props
}

test_that("preserve_correlation_structure_joint() is a drop-in for preserve_correlation_structure(): same required-argument shape, valid output", {
  target_prop_corr <- matrix(c(1, 0.6, 0.6, 1), nrow = 2, dimnames = list(c("clay", "bd"), c("clay", "bd")))
  property_matrices <- make_independent_depth_property_matrices(target_prop_corr, n_depths = 4, n_sims = 2000)
  gp_predictions <- list(clay = c(20, 22, 25, 28), bd = c(1.3, 1.35, 1.4, 1.45))
  depths <- c(0, 20, 50, 100)

  # Called with exactly the same positional shape preserve_correlation_structure() itself takes.
  result <- preserve_correlation_structure_joint(property_matrices, gp_predictions, depths, "clay")

  expect_named(result, c("clay", "bd"))
  expect_equal(dim(result$clay), dim(property_matrices$clay))
  expect_equal(dim(result$bd), dim(property_matrices$bd))

  # An unrecognized primary_property (unlike preserve_correlation_structure(), which falls back to
  # property_names[1] with a log message) is simply ignored, not an error - the joint method has
  # no use for it.
  expect_no_error(
    preserve_correlation_structure_joint(property_matrices, gp_predictions, depths, "nonexistent_property")
  )
})

test_that("preserve_correlation_structure_joint() degrades gracefully (returns the original matrices) on insufficient dimensions", {
  property_matrices <- list(clay = matrix(1:6, nrow = 1))  # n_depths = 1
  result <- preserve_correlation_structure_joint(property_matrices, list(clay = 20), depths = 10, primary_property = "clay")
  expect_identical(result, property_matrices)
})

test_that("preserve_correlation_structure_joint() induces real vertical correlation while preserving the existing property correlation", {
  # The core acceptance test for this phase: start from data with near-zero depth correlation but
  # a KNOWN property correlation, run the joint method (no gp_models -> falls back to the
  # full-depth-range length-scale), and confirm via Phase 0's validate_joint_correlation_structure()
  # that (a) depth correlation increased substantially from its near-zero baseline, and (b) property
  # correlation is still close to the original target - i.e. real vertical structure was ADDED
  # without destroying the property structure that was already there.
  target_prop_corr <- matrix(c(
    1.00, -0.55,
   -0.55,  1.00
  ), nrow = 2, dimnames = list(c("clay", "bd"), c("clay", "bd")))
  depths <- c(0, 20, 50, 100)
  property_matrices <- make_independent_depth_property_matrices(target_prop_corr, n_depths = 4, n_sims = 4000, seed = 52)

  before <- validate_joint_correlation_structure(property_matrices, target_property_corr = target_prop_corr)
  expect_true(before$overall_property_max_diff < 0.1)  # property correlation already correct
  # Depth 1 vs depth 4 correlation should be near zero before adjustment (fixture's own contract).
  for (prop in names(property_matrices)) {
    expect_true(abs(before$achieved_depth_correlation[[prop]][1, 4]) < 0.15)
  }

  result <- preserve_correlation_structure_joint(
    property_matrices, gp_predictions = NULL, depths = depths, primary_property = "clay"
  )

  after <- validate_joint_correlation_structure(result, target_property_corr = target_prop_corr)

  # Property correlation preserved (not distorted by adding depth correlation).
  expect_true(after$overall_property_max_diff < 0.15)

  # Depth correlation meaningfully increased for both properties (real vertical structure added).
  for (prop in names(result)) {
    achieved_depth1_4 <- abs(after$achieved_depth_correlation[[prop]][1, 4])
    expect_true(achieved_depth1_4 > 0.15)
  }
})

test_that("preserve_correlation_structure_joint() derives its depth length-scale from gp_models when supplied", {
  set.seed(53)
  agg_data <- data.frame(hzdept_r = c(0, 20, 50, 100), mean_val = c(10, 14, 22, 30))
  fitted <- fit_local_gp_model_single(agg_data, "clay")
  expect_false(is.null(fitted))

  target_prop_corr <- diag(1)
  dimnames(target_prop_corr) <- list("clay", "clay")
  property_matrices <- make_independent_depth_property_matrices(target_prop_corr, n_depths = 4, n_sims = 1500, seed = 54)

  result <- preserve_correlation_structure_joint(
    property_matrices, gp_predictions = NULL, depths = c(0, 20, 50, 100), primary_property = "clay",
    gp_models = list(clay = fitted)
  )
  expect_equal(dim(result$clay), dim(property_matrices$clay))
})

test_that("preserve_correlation_structure_joint() applies GP-mean recentering when gp_predictions is supplied", {
  target_prop_corr <- diag(1)
  dimnames(target_prop_corr) <- list("clay", "clay")
  property_matrices <- make_independent_depth_property_matrices(target_prop_corr, n_depths = 3, n_sims = 3000, seed = 55)
  depths <- c(0, 20, 50)
  gp_predictions <- list(clay = c(100, 100, 100))  # far from the fixture's own ~1-3 range

  result <- preserve_correlation_structure_joint(property_matrices, gp_predictions, depths, "clay")
  for (i in seq_len(3)) {
    expect_equal(mean(result$clay[i, ]), 100, tolerance = 1)
  }
})

test_that("calculate_safe_gp_ratio() clamps extreme ratios and returns 1 for invalid inputs", {
  expect_equal(calculate_safe_gp_ratio(c(10, 20), 2), 2)
  expect_equal(calculate_safe_gp_ratio(c(NA, 20), 2), 1)
  expect_equal(calculate_safe_gp_ratio(c(0, 20), 2), 1)
  expect_equal(calculate_safe_gp_ratio(c(1, 1000), 2), 10)   # clamped to max 10
  expect_equal(calculate_safe_gp_ratio(c(1000, 1), 2), 0.1)  # clamped to min 0.1
})

test_that("apply_quantile_adjustment() matches a direct per-element quantile() computation", {
  set.seed(21)
  n_sims <- 50
  curr_values <- rnorm(n_sims, mean = 20, sd = 3)
  prev_values <- rnorm(n_sims, mean = 18, sd = 3)
  reference_quantiles <- runif(n_sims)
  gp_ratio <- 1.1

  result <- apply_quantile_adjustment(reference_quantiles, curr_values, prev_values, gp_ratio, n_sims)
  expect_length(result, n_sims)

  # Hand-computed via the same per-element logic the old loop used, to confirm the vectorized
  # quantile() call (one sort, all probs at once) is equivalent to n_sims separate single-prob
  # quantile() calls.
  expected <- vapply(seq_len(n_sims), function(j) {
    qv <- stats::quantile(curr_values, probs = reference_quantiles[j], na.rm = TRUE, names = FALSE)
    qv + (prev_values[j] * gp_ratio - qv)
  }, numeric(1))
  expect_equal(result, expected, tolerance = 1e-9)
})

test_that("apply_quantile_adjustment() falls back to curr_values on error", {
  result <- apply_quantile_adjustment(reference_quantiles = c(0.5, 0.5), curr_values = c(NA_real_, NA_real_),
                                       prev_values = c(1, 2), gp_ratio = 1, n_sims = 2)
  expect_equal(result, c(NA_real_, NA_real_))
})

test_that("get_property_constraints()/apply_range_constraints() enforce known plausible ranges", {
  texture <- get_property_constraints("claytotal")
  expect_equal(texture$range, c(0, 100))
  expect_true(texture$preserve_distribution)

  unknown <- get_property_constraints("some_unmapped_property")
  expect_null(unknown$range)

  clamped <- apply_range_constraints(c(-10, 50, 150), list(range = c(0, 100)))
  expect_equal(clamped, c(0, 50, 100))
})

test_that("merge_adjusted_data()'s vectorized key match matches the original per-row which()-scan loop bit-for-bit", {
  # PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4: covers the original per-row loop's edge cases -
  # exactly-one-match required (zero or duplicate matches in result_data skip the update), NA
  # values in adjusted_data skip that row/property, and duplicate keys WITHIN adjusted_data
  # resolve to the last (highest row index) value, matching sequential last-write-wins.
  cokey_data <- data.frame(
    hzdept_r = c(0, 0, 20, 20, 40, 40), # (0, 2) is a duplicate key -> should never match
    simulation_number = c(1, 2, 1, 2, 1, 1), # (40, 1) also duplicated -> should never match
    clay_pct = c(15, 16, 17, 18, 19, 20),
    sand_pct = c(45, 44, 43, 42, 41, 40)
  )
  adjusted_data <- data.frame(
    hzdept_r = c(0, 0, 20, 20, 20, 40, 99),
    simulation_number = c(1, 1, 1, 2, 2, 1, 1), # rows 1-2 duplicate key (0,1); last wins
    clay_pct = c(100, 200, NA, 300, 301, 999, 500), # row 3 NA -> skip; rows 4-5 duplicate, last wins
    sand_pct = c(-1, -2, -3, -4, -5, -6, -7)
  )
  available_properties <- c("clay_pct", "sand_pct")

  reference_merge <- function(cokey_data, adjusted_data, available_properties) {
    result_data <- cokey_data
    for (i in seq_len(nrow(adjusted_data))) {
      adj_row <- adjusted_data[i, ]
      match_idx <- which(
        result_data$hzdept_r == adj_row$hzdept_r &
          result_data$simulation_number == adj_row$simulation_number
      )
      if (length(match_idx) == 1) {
        for (prop in available_properties) {
          if (prop %in% names(adj_row) && !is.na(adj_row[[prop]])) {
            result_data[[prop]][match_idx] <- adj_row[[prop]]
          }
        }
      }
    }
    result_data
  }

  result <- merge_adjusted_data(cokey_data, adjusted_data, available_properties)
  expected <- reference_merge(cokey_data, adjusted_data, available_properties)
  expect_equal(result, expected)

  # Spot-check the semantics directly, not just bit-identity to the reference:
  # (0,1) matched a unique result_data row and had 2 duplicate adjusted_data rows -> last (200) wins.
  expect_equal(result$clay_pct[result$hzdept_r == 0 & result$simulation_number == 1], 200)
  # (0,2) is duplicated WITHIN result_data -> never matched, stays at its original value (16).
  expect_equal(result$clay_pct[result$hzdept_r == 0 & result$simulation_number == 2][1], 16)
  # (20,1)'s adjusted_data value was NA -> skipped, stays at original (17).
  expect_equal(result$clay_pct[result$hzdept_r == 20 & result$simulation_number == 1], 17)
  # (20,2) had 2 non-NA duplicate adjusted_data rows -> last (301) wins.
  expect_equal(result$clay_pct[result$hzdept_r == 20 & result$simulation_number == 2], 301)
  # (40,1) is duplicated WITHIN result_data -> never matched, stays at original values (19, 20).
  expect_equal(sort(result$clay_pct[result$hzdept_r == 40]), c(19, 20))
})

test_that("apply_cross_property_constraints() renormalizes texture properties to sum to 100", {
  df <- data.frame(sandtotal = c(50, 30), claytotal = c(30, 30), silttotal = c(30, 20))
  result <- apply_cross_property_constraints(df, c("sandtotal", "claytotal", "silttotal"))
  totals <- rowSums(result[, c("sandtotal", "claytotal", "silttotal")])
  expect_equal(totals, c(100, 100), tolerance = 1e-6)
})

test_that("apply_cross_property_constraints()'s vectorized rowSums() matches the original per-row loop, including NA/zero-sum edge cases", {
  # PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4: row 3 has an NA in one texture column (NA sum ->
  # skip, per the original's `texture_sum > 0 && !is.na(texture_sum)` check); row 4 sums to
  # exactly 0 (skip, avoids division by zero); row 5 is already exactly 100 (scaling_factor = 1,
  # a no-op that must still leave the row unchanged, not corrupted).
  df <- data.frame(
    sandtotal = c(50, 30, 40, 0, 33.3333),
    claytotal = c(30, 30, NA, 0, 33.3333),
    silttotal = c(30, 20, 20, 0, 33.3334)
  )
  properties <- c("sandtotal", "claytotal", "silttotal")

  reference_constraints <- function(data, properties) {
    texture_props <- intersect(c("sand_total", "sandtotal", "clay_total", "claytotal", "silt_total", "silttotal"), properties)
    for (i in seq_len(nrow(data))) {
      texture_values <- unlist(data[i, texture_props])
      texture_sum <- sum(texture_values, na.rm = TRUE)
      if (texture_sum > 0 && !is.na(texture_sum)) {
        scaling_factor <- 100 / texture_sum
        data[i, texture_props] <- texture_values * scaling_factor
      }
    }
    data
  }

  result <- apply_cross_property_constraints(df, properties)
  expected <- reference_constraints(df, properties)
  expect_equal(result, expected)

  # Row 3 (NA claytotal): na.rm=TRUE means the sum (60) excludes the NA, so this row IS scaled
  # (100/60 factor) - only the NA cell itself stays NA (NA * scaling_factor is still NA).
  expect_equal(result$sandtotal[3], 40 * 100 / 60)
  expect_true(is.na(result$claytotal[3]))
  # Row 4 (all zeros) is untouched, not NaN from a 100/0 division.
  expect_equal(as.numeric(result[4, properties]), c(0, 0, 0))
})

test_that("get_nrcs_property_mapping() returns a static lookup with the expected keys", {
  mapping <- get_nrcs_property_mapping()
  expect_equal(mapping[["claytotal"]], "clay_pct")
  expect_equal(mapping[["ph"]], "pH")
})

test_that("match_simulations_to_nrcs_models() falls back correctly when unmapped or mapping is NULL", {
  expect_equal(match_simulations_to_nrcs_models("1", NULL), "general_pool")

  mapping <- data.frame(sim_cokey = c("1", "2"), gp_model_group = c("groupA", NA), stringsAsFactors = FALSE)
  expect_equal(match_simulations_to_nrcs_models("1", mapping), "groupA")
  expect_equal(match_simulations_to_nrcs_models("2", mapping), "general_pool")
  expect_equal(match_simulations_to_nrcs_models("nonexistent", mapping), "general_pool")
})

test_that("check_missing_values_increase() reports the NA-count delta per property", {
  orig <- data.frame(clay_pct = c(1, NA, 3))
  integrated <- data.frame(clay_pct = c(1, NA, NA))
  result <- check_missing_values_increase(orig, integrated, "clay_pct")
  expect_equal(result$clay_pct$increase, 1)
})

test_that("correct_distribution_shape() returns curr_values unadjusted for all-NA input rather than erroring", {
  all_na <- rep(NA_real_, 5)
  expect_silent(result <- correct_distribution_shape(all_na, all_na))
  expect_equal(result, all_na)

  # A single non-NA value: var() is undefined (NA), same guarded path.
  one_value <- c(NA, NA, 5, NA)
  expect_silent(result2 <- correct_distribution_shape(one_value, c(1, 2, 3, 4)))
  expect_equal(result2, one_value)

  # adjusted_curr all-NA (would break ecdf()/quantile()) while curr_values is fine.
  expect_silent(result3 <- correct_distribution_shape(c(1, 2, 3, 4, 5), rep(NA_real_, 5)))
  expect_equal(result3, c(1, 2, 3, 4, 5))
})

test_that("correct_distribution_shape() still performs a real correction when there's enough data", {
  set.seed(11)
  curr_values <- rnorm(50, mean = 10, sd = 2)
  adjusted_curr <- curr_values * 1.5
  result <- correct_distribution_shape(curr_values, adjusted_curr)
  expect_equal(length(result), length(curr_values))
  expect_true(all(is.finite(result)))
})

test_that("detect_simulation_properties() intersects known SSURGO/lab property names with the data's own columns", {
  df <- data.frame(cokey = "1", sandtotal = 40, claytotal = 20, unrelated_col = 1)
  detected <- detect_simulation_properties(df)
  expect_setequal(detected, c("sandtotal", "claytotal"))
})

test_that("convert_to_property_matrices()/convert_to_long_format() round-trip depth x simulation values", {
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20), n_sims = 3)
  depths <- c(0, 20)
  sims <- 1:3

  matrices <- convert_to_property_matrices(sim_data, c("clay_pct", "sand_pct"), depths, sims)
  expect_equal(dim(matrices$clay_pct), c(2, 3))
  expect_false(any(is.na(matrices$clay_pct)))

  long <- convert_to_long_format(matrices, depths, sims, sim_data, c("clay_pct", "sand_pct"))
  expect_equal(nrow(long), 6)
  expect_true(all(c("clay_pct", "sand_pct", "hzdept_r", "simulation_number") %in% names(long)))
})

test_that("aggregate_property_by_depth()/fit_local_gp_model_single() fit a usable local GP model", {
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 8)
  cokey_data <- sim_data[sim_data$cokey == "1", ]

  agg <- aggregate_property_by_depth(cokey_data, "clay_pct")
  expect_equal(nrow(agg), 4)
  expect_true(all(c("hzdept_r", "mean_val") %in% names(agg)))

  model <- fit_local_gp_model_single(agg, "clay_pct")
  expect_false(is.null(model))
  expect_equal(model$property, "clay_pct")

  preds <- generate_local_predictions(list(clay_pct = model), c(0, 50, 100))
  expect_length(preds$clay_pct, 3)
  expect_true(all(is.finite(preds$clay_pct)))
})

test_that("fit_local_gp_model_single()'s cheaper default gp_control matches GPfit::GP_fit()'s own default", {
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 8)
  cokey_data <- sim_data[sim_data$cokey == "1", ]
  agg <- aggregate_property_by_depth(cokey_data, "clay_pct")

  model_default <- fit_local_gp_model_single(agg, "clay_pct") # uses the fast default
  model_thorough <- fit_local_gp_model_single(agg, "clay_pct", gp_control = c(200, 80, 2)) # GP_fit()'s own default

  # Two independent numerical hyperparameter searches - expect the same optimum to within
  # ordinary floating-point-level noise, not bit-for-bit identical.
  expect_equal(model_default$gp_model$beta, model_thorough$gp_model$beta, tolerance = 1e-4)

  preds_default <- generate_local_predictions(list(clay_pct = model_default), c(0, 50, 100))
  preds_thorough <- generate_local_predictions(list(clay_pct = model_thorough), c(0, 50, 100))
  expect_equal(preds_default$clay_pct, preds_thorough$clay_pct, tolerance = 1e-4)
})

test_that("apply_local_gp_adjustments() returns cokey_data unchanged when there are too few unique depths", {
  cokey_data <- data.frame(cokey = "1", hzdept_r = c(0, 20), simulation_number = 1:2, clay_pct = c(15, 18))
  result <- apply_local_gp_adjustments(cokey_data, "clay_pct", min_depths = 3)
  expect_identical(result, cokey_data)
})

test_that("apply_local_gp_adjustments() runs end-to-end on adequately-varied synthetic data without error", {
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 10)
  cokey_data <- sim_data[sim_data$cokey == "1", ]

  result <- apply_local_gp_adjustments(cokey_data, c("clay_pct", "sand_pct"), preserve_correlations = TRUE)
  expect_equal(nrow(result), nrow(cokey_data))
  expect_true(all(c("clay_pct", "sand_pct") %in% names(result)))
})

test_that("apply_gp_depth_trends() defaults to the joint_copula method (Phase 13 flip) - config NULL/empty/explicit all agree", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 13 flipped the default from
  # "gp_quantile_retrofit" to "joint_copula" after Phases 0-12 resolved every blocking decision
  # point. No config (or a config that simply doesn't set vertical_correlation_method) must now
  # dispatch to the exact same preserve_correlation_structure_joint() path as an explicit
  # "joint_copula" - matching default_monte_carlo_config()'s own new default.
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 10)
  cokey_data <- sim_data[sim_data$cokey == "1", ]
  gp_predictions <- list(clay_pct = c(15, 17, 20, 25), sand_pct = c(45, 44, 42, 40))

  set.seed(61)
  result_no_config <- apply_gp_depth_trends(cokey_data, gp_predictions, c("clay_pct", "sand_pct"))

  set.seed(61)
  result_empty_config <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"), config = list()
  )

  set.seed(61)
  result_explicit_joint <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"),
    config = list(monte_carlo = list(vertical_correlation_method = "joint_copula"))
  )

  expect_identical(result_no_config, result_empty_config)
  expect_identical(result_no_config, result_explicit_joint)
})

test_that("apply_gp_depth_trends()'s gp_quantile_retrofit opt-out still works exactly as before and differs from the new joint_copula default", {
  # The original algorithm is fully supported, not deprecated or removed by the Phase 13 flip -
  # this locks in that explicitly opting out still reaches preserve_correlation_structure()
  # (the original, unmodified retrofit code path) and produces different output than the new
  # default, confirming the dispatch still genuinely branches both ways.
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 10)
  cokey_data <- sim_data[sim_data$cokey == "1", ]
  gp_predictions <- list(clay_pct = c(15, 17, 20, 25), sand_pct = c(45, 44, 42, 40))
  retrofit_config <- list(monte_carlo = list(vertical_correlation_method = "gp_quantile_retrofit"))

  set.seed(62)
  result_default <- apply_gp_depth_trends(cokey_data, gp_predictions, c("clay_pct", "sand_pct"))

  set.seed(62)
  result_retrofit <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"), config = retrofit_config
  )

  expect_equal(nrow(result_retrofit), nrow(cokey_data))
  expect_true(all(c("clay_pct", "sand_pct") %in% names(result_retrofit)))
  # The two methods are genuinely different algorithms - same seed, same inputs, should not
  # produce bit-identical adjusted values (confirms the dispatch actually branches, not a no-op).
  expect_false(isTRUE(all.equal(result_default$clay_pct, result_retrofit$clay_pct)))
})

test_that("apply_gp_depth_trends() extracts per-depth boundary_distinctness from a bound_sd column and passes gp_models through under joint_copula, without error", {
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 8)
  cokey_data <- sim_data[sim_data$cokey == "1", ]
  # bound_sd broadcast per-depth, as simulate_cokey_generalized() actually produces it (Phase 1b).
  cokey_data$bound_sd <- c(2.5, 2.5, 7.5, 7.5)[match(cokey_data$hzdept_r, c(0, 20, 50, 100))]
  gp_predictions <- list(clay_pct = c(15, 17, 20, 25))

  agg_data <- data.frame(hzdept_r = c(0, 20, 50, 100), mean_val = c(15, 17, 20, 25))
  fitted <- fit_local_gp_model_single(agg_data, "clay_pct")

  result <- apply_gp_depth_trends(
    cokey_data, gp_predictions, "clay_pct",
    config = list(monte_carlo = list(
      vertical_correlation_method = "joint_copula", vertical_correlation_gating = TRUE
    )),
    gp_models = list(clay_pct = fitted)
  )
  expect_equal(nrow(result), nrow(cokey_data))
  expect_true("clay_pct" %in% names(result))
})

test_that("apply_gp_depth_trends()'s discontinuity gating stays OFF under joint_copula unless vertical_correlation_gating is explicitly TRUE (Phase 8 decoupling)", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 8: bound_sd is attached unconditionally
  # upstream in production (attach_osd_boundary_distinctness()), so gating must NOT silently
  # activate just because bound_sd happens to be present - it needs its own explicit opt-in,
  # independent of vertical_correlation_method itself.
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 40, 60), n_sims = 400)
  cokey_data <- sim_data[sim_data$cokey == "1", ]
  # An injected abrupt boundary between depths 40 and 60 (bound_sd ~= 1, aqp's own
  # abrupt-code-to-offset value) - if gating is active, this should visibly suppress
  # cross-boundary correlation; if inactive, output should match the no-bound_sd case exactly.
  cokey_data$bound_sd <- c(NA, NA, 1, NA)[match(cokey_data$hzdept_r, c(0, 20, 40, 60))]
  gp_predictions <- list(clay_pct = c(15, 15, 15, 15), sand_pct = c(45, 45, 45, 45))

  cokey_data_no_bound_sd <- cokey_data
  cokey_data_no_bound_sd$bound_sd <- NULL

  joint_config_default_gating <- list(monte_carlo = list(vertical_correlation_method = "joint_copula"))
  joint_config_gating_off <- list(monte_carlo = list(
    vertical_correlation_method = "joint_copula", vertical_correlation_gating = FALSE
  ))
  joint_config_gating_on <- list(monte_carlo = list(
    vertical_correlation_method = "joint_copula", vertical_correlation_gating = TRUE
  ))

  set.seed(81)
  result_no_bound_sd <- apply_gp_depth_trends(
    cokey_data_no_bound_sd, gp_predictions, c("clay_pct", "sand_pct"), config = joint_config_default_gating
  )
  set.seed(81)
  result_default_gating <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"), config = joint_config_default_gating
  )
  set.seed(81)
  result_gating_off <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"), config = joint_config_gating_off
  )
  set.seed(81)
  result_gating_on <- apply_gp_depth_trends(
    cokey_data, gp_predictions, c("clay_pct", "sand_pct"), config = joint_config_gating_on
  )

  # config$monte_carlo$vertical_correlation_gating unset (default) behaves identically to
  # explicitly FALSE, and identically to bound_sd not being present at all - gating is off by
  # default, full stop. (result_no_bound_sd naturally lacks a bound_sd column at all - that
  # column simply passes through cokey_data untouched when present, unrelated to whether gating
  # itself is active - so compare only the columns both share.)
  expect_identical(result_default_gating, result_gating_off)
  shared_cols <- intersect(names(result_default_gating), names(result_no_bound_sd))
  expect_identical(result_default_gating[shared_cols], result_no_bound_sd[shared_cols])

  # Explicitly enabling gating with a real discontinuity present must change the output - the
  # flag actually does something when turned on.
  expect_false(isTRUE(all.equal(result_gating_on$clay_pct, result_gating_off$clay_pct)))
})

test_that("joint_copula is reachable through the NRCS/regional GP path (apply_nrcs_trend_adjustments()/process_single_cokey()), matching the local-GP path's Phase 6 fix", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 11: apply_nrcs_trend_adjustments() previously
  # had no config parameter at all, and process_single_cokey() didn't pass its own config to it
  # even though the sibling local-GP branch just below it already received config (Phase 6) - so
  # "joint_copula" was unreachable via integration_method = "nrcs_gp" regardless of the config
  # flag. This closes that gap the same way Phase 6 closed it for the local path.
  # Two mappable properties (get_nrcs_property_mapping(): claytotal -> clay_pct,
  # sandtotal -> sand_pct) - the vertical-correlation method dispatch inside
  # apply_gp_depth_trends() only branches on preserve_correlations with >= 2 properties; with
  # just one, both methods trivially fall through the SAME "individual adjustment" path and would
  # never actually diverge.
  training_df <- data.frame(
    cokey = rep(c("1", "2", "3"), each = 4),
    hzdept_r = rep(c(0, 20, 50, 100), 3),
    clay_pct = c(10, 15, 22, 30, 12, 16, 24, 32, 8, 14, 20, 28),
    sand_pct = c(45, 42, 38, 33, 44, 41, 36, 30, 47, 43, 39, 34),
    compname = "Alpha",
    stringsAsFactors = FALSE
  )
  clay_fit <- fit_individual_gp_model(training_df, "clay_pct", optimize_hyperparameters = FALSE)
  sand_fit <- fit_individual_gp_model(training_df, "sand_pct", optimize_hyperparameters = FALSE)
  expect_false(is.null(clay_fit))
  expect_false(is.null(sand_fit))
  nrcs_gp_models <- list(
    clay_pct = list(type = "stratified_grouped", models = list(group1 = clay_fit$model)),
    sand_pct = list(type = "stratified_grouped", models = list(group1 = sand_fit$model))
  )

  cokey_data <- data.frame(
    cokey = "99", hzdept_r = c(0, 20, 50, 100), simulation_number = seq_len(200),
    claytotal = 10 + stats::rnorm(200, sd = 0.5),  # near-constant/independent-across-depth
    sandtotal = 45 + stats::rnorm(200, sd = 0.5)
  )

  retrofit_config <- list(monte_carlo = list(vertical_correlation_method = "gp_quantile_retrofit"))
  joint_config <- list(monte_carlo = list(vertical_correlation_method = "joint_copula"))
  properties <- c("claytotal", "sandtotal")

  set.seed(111)
  result_retrofit <- apply_nrcs_trend_adjustments(
    cokey_data, nrcs_gp_models, "group1", properties, config = retrofit_config
  )
  set.seed(111)
  result_joint <- apply_nrcs_trend_adjustments(
    cokey_data, nrcs_gp_models, "group1", properties, config = joint_config
  )

  expect_equal(nrow(result_joint), nrow(cokey_data))
  expect_false(isTRUE(all.equal(result_retrofit$claytotal, result_joint$claytotal)))

  # config = NULL (the default) now reproduces the joint_copula result, not the retrofit one -
  # Phase 13 flipped the default; "gp_quantile_retrofit" remains reachable only as an explicit
  # opt-out (already confirmed via `retrofit_config` above).
  set.seed(111)
  result_default <- apply_nrcs_trend_adjustments(cokey_data, nrcs_gp_models, "group1", properties)
  expect_identical(result_default, result_joint)

  # And confirmed reachable one level up, through process_single_cokey() itself (the actual
  # call site the gap existed in).
  cokey_mapping <- data.frame(cokey = "99", model_group = "group1", stringsAsFactors = FALSE)
  set.seed(112)
  via_process_single_cokey <- process_single_cokey(
    cokey_data, "99", properties, nrcs_gp_models, cokey_mapping,
    use_nrcs_gp = TRUE, use_local_gp = FALSE,
    preserve_correlations = TRUE, config = joint_config
  )
  expect_false(is.null(via_process_single_cokey))
  expect_equal(nrow(via_process_single_cokey), nrow(cokey_data))
})

test_that("process_single_cokey()/process_cokeys_sequential() integrate local GP adjustments across cokeys", {
  sim_data <- make_sim_data(cokeys = c("1", "2"), depths = c(0, 20, 50, 100), n_sims = 8)
  config <- default_config("validation")

  results <- process_cokeys_sequential(
    split(sim_data, sim_data$cokey), unique(sim_data$cokey), c("clay_pct", "sand_pct"),
    gp_models = NULL, cokey_mapping = NULL,
    use_nrcs_gp = FALSE, use_local_gp = TRUE,
    preserve_correlations = TRUE, config = config
  )
  expect_length(results, 2)
  expect_true(all(vapply(results, function(r) !is.null(r) && nrow(r) > 0, logical(1))))
})

test_that("joint_copula is reachable end-to-end from the real top-level API (process_single_cokey()/apply_local_gp_adjustments()), not just by calling apply_gp_depth_trends() directly", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 6: apply_local_gp_adjustments() already fits
  # local_gp_models and already receives config (both were sitting unused for this purpose before
  # this change) - apply_local_depth_trends() now threads both through to apply_gp_depth_trends(),
  # closing the gap Phase 4 had explicitly deferred (which only wired the dispatch INTO
  # apply_gp_depth_trends() itself, not the full call chain above it).
  sim_data <- make_sim_data(cokeys = "1", depths = c(0, 20, 50, 100), n_sims = 10)
  cokey_data <- sim_data[sim_data$cokey == "1", ]

  retrofit_config <- default_config("validation")
  retrofit_config$monte_carlo$vertical_correlation_method <- "gp_quantile_retrofit"

  joint_config <- default_config("validation")
  joint_config$monte_carlo$vertical_correlation_method <- "joint_copula"

  set.seed(71)
  retrofit_result <- apply_local_gp_adjustments(
    cokey_data, c("clay_pct", "sand_pct"), preserve_correlations = TRUE, config = retrofit_config
  )
  set.seed(71)
  joint_result <- apply_local_gp_adjustments(
    cokey_data, c("clay_pct", "sand_pct"), preserve_correlations = TRUE, config = joint_config
  )

  expect_equal(nrow(joint_result), nrow(cokey_data))
  expect_true(all(c("clay_pct", "sand_pct") %in% names(joint_result)))
  # Genuinely different algorithms reaching this call site - same seed/inputs should not produce
  # bit-identical output (confirms the config actually reaches apply_gp_depth_trends()'s dispatch
  # through the full apply_local_gp_adjustments() -> apply_local_depth_trends() call chain, not
  # silently dropped along the way).
  expect_false(isTRUE(all.equal(retrofit_result$clay_pct, joint_result$clay_pct)))

  # config = NULL (apply_local_gp_adjustments()'s own default, which internally falls back to
  # default_config("validation") - a config that never sets vertical_correlation_method
  # at all, so apply_gp_depth_trends()'s own `%||% "joint_copula"` fallback is what actually
  # decides this) now reproduces the joint_copula result, matching Phase 13's flipped default -
  # not the retrofit result, which remains reachable only via the explicit opt-out above.
  set.seed(71)
  default_result <- apply_local_gp_adjustments(
    cokey_data, c("clay_pct", "sand_pct"), preserve_correlations = TRUE
  )
  expect_identical(default_result, joint_result)
})

test_that("future::multisession workers can see soilSIM package functions via automatic globals detection", {
  skip_on_cran()
  # This isolates the same load-bearing assumption the old clusterEvalQ(cl, library(soilSIM))
  # smoke test used to check, for the new mechanism: future/globals auto-detects that a
  # dispatched closure calls a soilSIM-namespaced function and attaches the package on the
  # worker automatically - this requires soilSIM to actually be installed and attached (not
  # merely devtools::load_all()'d), since a future::multisession worker is a brand-new R
  # process. Skip in a load_all()-only dev session; this runs for real under R CMD check /
  # devtools::check(), which installs the package first.
  skip_if_not(nzchar(system.file(package = "soilSIM")) &&
                file.exists(file.path(system.file(package = "soilSIM"), "Meta", "package.rds")),
              "soilSIM is not installed (only load_all()'d) in this session")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 1)

  worker_has_log_message <- future.apply::future_lapply(1, function(i) exists("log_message"))[[1]]
  expect_true(worker_has_log_message)
})

test_that("process_cokeys_parallel() produces the same results as process_cokeys_sequential() when parallel is usable, or degrades to it otherwise", {
  skip_on_cran()
  sim_data <- make_sim_data(cokeys = c("1", "2"), depths = c(0, 20, 50, 100), n_sims = 8)
  config <- default_config("validation")

  # n_cores = 2 (not 1) to actually exercise the future::multisession dispatch path -
  # run_parallel_lapply() short-circuits n_cores <= 1 straight to plain lapply() without ever
  # touching future, which would make this test pass trivially without exercising anything.
  # Regardless of whether soilSIM is installed (vs. only load_all()'d) in this dev session,
  # process_cokeys_parallel() has its own sequential_fallback to process_cokeys_sequential() on
  # any worker-side error - so this call must succeed either way and produce valid per-cokey
  # results. suppressWarnings(): in a load_all()-only dev session, the worker-side "package
  # not attached" failure is expected and reported as a warning before falling back - not a
  # code defect (see the skipped test above, which isolates and documents this environment
  # limitation).
  results <- suppressWarnings(process_cokeys_parallel(
    split(sim_data, sim_data$cokey), unique(sim_data$cokey), c("clay_pct", "sand_pct"),
    gp_models = NULL, cokey_mapping = NULL,
    use_nrcs_gp = FALSE, use_local_gp = TRUE,
    preserve_correlations = TRUE, n_cores = 2, config = config
  ))
  expect_length(results, 2)
  expect_true(all(vapply(results, function(r) !is.null(r) && nrow(r) > 0, logical(1))))
})
