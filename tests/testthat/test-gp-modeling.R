## Synthetic GP-fit-shaped data throughout - GPfit::GP_fit()/predict.GP() are
## pure local computation, so no network mocking is needed here.

make_gp_training_df <- function() {
  data.frame(
    cokey = rep(c("1", "2", "3"), each = 4),
    hzdept_r = rep(c(0, 20, 50, 100), 3),
    clay_pct = c(10, 15, 22, 30,
                 12, 16, 24, 32,
                 8, 14, 20, 28),
    compname = "Alpha",
    stringsAsFactors = FALSE
  )
}

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

test_that("adjust_multivariate_depthwise_GP() preserves matrix dimensions and validates mismatched inputs", {
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
  # adjust_multivariate_depthwise_GP() calls predict_gp_depth_trends() when the
  # element has a $gp_model - stub it here via a minimal fake to avoid needing
  # a real fitted GPfit object for this dimension/validation-focused test.
  testthat::local_mocked_bindings(
    predict_gp_depth_trends = function(gp_model_info, new_depths) gp_model_info$predictions_override,
    .package = "soilSIM"
  )

  result <- adjust_multivariate_depthwise_GP(simulated_list, gp_models, depths)
  expect_named(result, c("clay_pct", "sand_pct"))
  expect_equal(dim(result$clay_pct), c(n_depths, n_sims))
  expect_equal(dim(result$sand_pct), c(n_depths, n_sims))

  expect_error(adjust_multivariate_depthwise_GP(list(), gp_models, depths), "cannot be empty")
  expect_error(adjust_multivariate_depthwise_GP(simulated_list, gp_models, c(0, 50)), "Length of depths")
})

test_that("adjust_multivariate_depthwise_GP()'s vectorized quantile() call matches the original per-replicate loop", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 3 adjust_multivariate_depthwise_GP()
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
  actual <- adjust_multivariate_depthwise_GP(simulated_list, gp_models, depths, primary_property = "clay_pct")

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
  # adjust_multivariate_depthwise_GP() itself warns on.
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
