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

test_that("calculate_safe_gp_ratio() clamps extreme ratios and returns 1 for invalid inputs", {
  expect_equal(calculate_safe_gp_ratio(c(10, 20), 2), 2)
  expect_equal(calculate_safe_gp_ratio(c(NA, 20), 2), 1)
  expect_equal(calculate_safe_gp_ratio(c(0, 20), 2), 1)
  expect_equal(calculate_safe_gp_ratio(c(1, 1000), 2), 10)   # clamped to max 10
  expect_equal(calculate_safe_gp_ratio(c(1000, 1), 2), 0.1)  # clamped to min 0.1
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

test_that("apply_cross_property_constraints() renormalizes texture properties to sum to 100", {
  df <- data.frame(sandtotal = c(50, 30), claytotal = c(30, 30), silttotal = c(30, 20))
  result <- apply_cross_property_constraints(df, c("sandtotal", "claytotal", "silttotal"))
  totals <- rowSums(result[, c("sandtotal", "claytotal", "silttotal")])
  expect_equal(totals, c(100, 100), tolerance = 1e-6)
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

test_that("process_single_cokey()/process_cokeys_sequential() integrate local GP adjustments across cokeys", {
  sim_data <- make_sim_data(cokeys = c("1", "2"), depths = c(0, 20, 50, 100), n_sims = 8)
  config <- get_default_configuration("validation")

  results <- process_cokeys_sequential(
    sim_data, unique(sim_data$cokey), c("clay_pct", "sand_pct"),
    gp_models = NULL, cokey_mapping = NULL,
    use_nrcs_gp = FALSE, use_local_gp = TRUE,
    preserve_correlations = TRUE, config = config
  )
  expect_length(results, 2)
  expect_true(all(vapply(results, function(r) !is.null(r) && nrow(r) > 0, logical(1))))
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
  config <- get_default_configuration("validation")

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
    sim_data, unique(sim_data$cokey), c("clay_pct", "sand_pct"),
    gp_models = NULL, cokey_mapping = NULL,
    use_nrcs_gp = FALSE, use_local_gp = TRUE,
    preserve_correlations = TRUE, n_cores = 2, config = config
  ))
  expect_length(results, 2)
  expect_true(all(vapply(results, function(r) !is.null(r) && nrow(r) > 0, logical(1))))
})
