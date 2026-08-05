## Synthetic workflow-shaped fixtures throughout - no network required.

test_that("assess_correlation_differences() computes real max/mean/RMSE differences between two matrices", {
  orig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  sim <- matrix(c(1, 0.3, 0.3, 1), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  result <- assess_correlation_differences(orig, sim, criteria = list())
  expect_equal(result$max_difference, 0.2, tolerance = 1e-8)
  expect_equal(result$mean_difference, 0.1, tolerance = 1e-8)
})

test_that("get_realistic_property_ranges()/get_default_property_constraints() return usable, consistent constraint tables", {
  ranges <- get_realistic_property_ranges()
  expect_equal(ranges$ph, list(min = 2.5, max = 11.0))
  expect_true(all(vapply(ranges, function(r) r$min < r$max, logical(1))))

  constraints <- get_default_property_constraints()
  expect_identical(constraints$property_ranges, ranges)
  expect_true(constraints$enhanced_checks$use_module8_validation)
})

test_that("determine_quality_grade() maps scores to the documented letter-grade bands", {
  expect_equal(determine_quality_grade(0.95), "Excellent")
  expect_equal(determine_quality_grade(0.85), "Good")
  expect_equal(determine_quality_grade(0.75), "Acceptable")
  expect_equal(determine_quality_grade(0.65), "Needs Improvement")
  expect_equal(determine_quality_grade(0.3), "Poor")
})

test_that("calculate_weighted_quality_score() weights component scores (numeric weight vector) and falls back to a neutral score when empty", {
  scores <- list(monte_carlo = 1.0, correlation = 0.5)
  # weighted.mean() requires a numeric weights vector - a named list (as
  # produced by get_default_quality_weights(), see the dedicated quirk test
  # below) does not work here.
  weights <- c(monte_carlo = 0.75, correlation = 0.25)
  result <- calculate_weighted_quality_score(scores, weights)
  expect_equal(result, 1.0 * 0.75 + 0.5 * 0.25, tolerance = 1e-8)

  expect_equal(calculate_weighted_quality_score(list(), weights), 0.5)
})

test_that("QUIRK: calculate_weighted_quality_score() silently falls back to a simple mean when weights is a list (as get_default_quality_weights() returns), not a numeric vector", {
  # Pre-existing behavior, not something to "fix" as part of this migration:
  # weighted.mean() rejects list weights with "non-numeric argument to binary
  # operator"; the function's own tryCatch() catches that and falls back to
  # an unweighted mean() - logged as an ERROR+WARN, not silently swallowed,
  # but the caller still gets a (differently-computed) numeric result rather
  # than a hard failure. Since assess_workflow_quality() calls this with
  # get_default_quality_weights()'s list output whenever validation_config
  # doesn't override quality_weights, this fallback path is in fact the
  # ordinary/current behavior, not a rare edge case.
  scores <- list(monte_carlo = 1.0, correlation = 0.5)
  weights_as_list <- get_default_quality_weights()
  result <- suppressWarnings(calculate_weighted_quality_score(scores, weights_as_list))
  expect_equal(result, mean(c(1.0, 0.5)))
})

test_that("calculate_confidence_level() rewards high, consistent scores and handles the empty case", {
  expect_equal(calculate_confidence_level(list()), 0.5)
  high_consistent <- calculate_confidence_level(list(a = 0.9, b = 0.9, c = 0.9))
  low_variable <- calculate_confidence_level(list(a = 0.9, b = 0.1, c = 0.5))
  expect_true(high_consistent > low_variable)
  expect_true(high_consistent <= 1 && high_consistent >= 0)
})

test_that("identify_critical_issues() flags validation failures and sub-threshold overall scores", {
  validation_results <- list(
    monte_carlo_validation = list(validation_failed = TRUE),
    correlation_validation = list(validation_failed = FALSE),
    overall_assessment = list(quality_score = 0.4)
  )
  issues <- identify_critical_issues(validation_results, list())
  expect_true(any(grepl("monte_carlo_validation", issues)))
  expect_true(any(grepl("critical threshold", issues)))
  expect_false(any(grepl("correlation_validation", issues)))
})

test_that("calculate_monte_carlo_score()/calculate_correlation_score() aggregate real sub-scores from nested validation results", {
  mc_score <- calculate_monte_carlo_score(list(
    convergence_assessment = list(converged = TRUE),
    coverage_assessment = list(overall_coverage = 90),
    distribution_fidelity = list(overall_fidelity = 0.8)
  ))
  expect_equal(mc_score, mean(c(1.0, 0.9, 0.8)), tolerance = 1e-8)

  corr_score <- calculate_correlation_score(list(
    preservation_assessment = list(overall_preservation = list(correlation_rmse = 0.05)),
    matrix_quality = list(overall_quality = 0.9)
  ))
  expect_equal(corr_score, mean(c(1 - 0.05 * 10, 0.9)), tolerance = 1e-8)
})

test_that("initialize_validation_structure()/extract_and_validate_components() build the expected scaffolding", {
  structure <- initialize_validation_structure(Sys.time(), list(), NULL)
  expect_true(all(c("monte_carlo_validation", "overall_assessment", "recommendations") %in% names(structure)))

  workflow_results <- list(
    simulation_data = data.frame(cokey = "1", simulation_number = 1),
    gp_models = list(clay_pct = list(type = "pooled")),
    correlation_matrices = list(global_correlation_matrix = diag(2))
  )
  components <- extract_and_validate_components(workflow_results)
  expect_true(is.data.frame(components$simulation_data))
  expect_false(is.null(components$gp_models))
  expect_false(is.null(components$correlation_matrices))
})

test_that("assess_workflow_quality() orchestrates real weighting/grading over component stub scores without erroring", {
  validation_results <- list(
    monte_carlo_validation = list(
      convergence_assessment = list(converged = TRUE),
      coverage_assessment = list(overall_coverage = 90),
      distribution_fidelity = list(overall_fidelity = 0.8)
    ),
    correlation_validation = list(
      preservation_assessment = list(overall_preservation = list(correlation_rmse = 0.02)),
      matrix_quality = list(overall_quality = 0.95)
    )
  )
  # suppressWarnings(): validation_config = list() means quality_weights
  # falls back to get_default_quality_weights()'s list output, which hits
  # the calculate_weighted_quality_score() list-weights quirk documented
  # above (logs ERROR+WARN, then falls back to an unweighted mean - not a
  # hard failure).
  assessment <- suppressWarnings(assess_workflow_quality(validation_results, list()))
  expect_true(assessment$quality_score >= 0 && assessment$quality_score <= 1)
  expect_true(assessment$quality_grade %in% c("Excellent", "Good", "Acceptable", "Needs Improvement", "Poor"))
  expect_true(assessment$workflow_status %in% c("PASSED", "FAILED"))
})

## ----------------------------------------------------------------------
## 4a. Monte Carlo / distribution assessment cluster
## ----------------------------------------------------------------------

test_that("assess_simulation_convergence() detects stabilized vs drifting realizations", {
  stable_data <- data.frame(simulation_number = 1:20, clay_pct = rep(20, 20) + stats::rnorm(20, 0, 0.01))
  result_stable <- assess_simulation_convergence(stable_data, list(tolerance = 0.05, n_batches = 4))
  expect_true(isTRUE(result_stable$converged))

  drifting_data <- data.frame(simulation_number = 1:20, clay_pct = seq(0, 100, length.out = 20))
  result_drifting <- assess_simulation_convergence(drifting_data, list(tolerance = 0.05, n_batches = 4))
  expect_false(isTRUE(result_drifting$converged))
})

test_that("assess_property_coverage() measures range overlap and flags extrapolation", {
  orig_values <- seq(10, 30, length.out = 50)
  good_sim <- seq(12, 28, length.out = 50)
  result_good <- assess_property_coverage(good_sim, orig_values, list())
  expect_true(result_good$coverage_percentage > 50)
  expect_false(result_good$extrapolation_detected)

  extrapolating_sim <- c(good_sim, 1000)
  result_bad <- assess_property_coverage(extrapolating_sim, orig_values, list(max_extrapolation_factor = 1.2))
  expect_true(result_bad$extrapolation_detected)
})

test_that("assess_distributional_coverage() flags high outlier percentages", {
  set.seed(1)
  clean_values <- stats::rnorm(200, 20, 2)
  result_clean <- assess_distributional_coverage(clean_values, list(max_outlier_percentage = 10))
  expect_true(result_clean$coverage_adequate)

  contaminated_values <- c(stats::rnorm(190, 20, 1), rep(1000, 10))
  result_contaminated <- assess_distributional_coverage(contaminated_values, list(max_outlier_percentage = 2))
  expect_false(result_contaminated$coverage_adequate)
})

test_that("assess_group_representation() flags underrepresented cokeys", {
  sim_data <- data.frame(cokey = c(rep("1", 15), rep("2", 3)))
  result <- assess_group_representation(sim_data, list(min_samples_per_group = 10))
  expect_equal(result$underrepresented_groups, "2")
  expect_false(result$representation_adequate)
})

test_that("test_distribution_fidelity()/compare_distribution_moments()/compare_distribution_quantiles() compare simulated values against a fitted normal distribution", {
  set.seed(1)
  param_info <- list(family = "normal", fit = list(mean = 20, sd = 2))
  matching_values <- stats::rnorm(300, 20, 2)

  fidelity <- test_distribution_fidelity(matching_values, param_info, list(ks_test_alpha = 0.05))
  expect_true(is.finite(fidelity$ks_test_p_value))

  moments <- compare_distribution_moments(matching_values, param_info, list(moment_tolerance = 0.2))
  expect_true(moments$mean_difference < 0.2)

  quantiles <- compare_distribution_quantiles(matching_values, param_info, list(quantile_tolerance = 0.2))
  expect_true(all(is.finite(quantiles$quantile_differences)))

  mismatched_values <- stats::rnorm(300, 100, 2)
  fidelity_mismatch <- test_distribution_fidelity(mismatched_values, param_info, list(ks_test_alpha = 0.05))
  expect_false(fidelity_mismatch$distribution_match)
})

## ----------------------------------------------------------------------
## 4b. Correlation-matrix quality cluster
## ----------------------------------------------------------------------

test_that("assess_correlation_matrix_properties() detects positive-definite vs non-positive-definite matrices", {
  pd_matrix <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  result_pd <- assess_correlation_matrix_properties(pd_matrix)
  expect_true(result_pd$is_positive_definite)
  expect_true(is.finite(result_pd$condition_number))

  non_pd_matrix <- matrix(c(1, 2, 2, 1), nrow = 2)  # eigenvalues -1, 3 - not PD
  result_bad <- assess_correlation_matrix_properties(non_pd_matrix)
  expect_false(result_bad$is_positive_definite)
})

test_that("validate_correlation_matrix_quality() aggregates across a named list of matrices", {
  matrices <- list(a = matrix(c(1, 0.2, 0.2, 1), nrow = 2), b = matrix(c(1, 2, 2, 1), nrow = 2))
  result <- validate_correlation_matrix_quality(matrices, list())
  expect_false(result$positive_definite)  # matrix "b" is not PD
})

test_that("assess_single_cholesky_decomposition() succeeds for a valid PD matrix and fails for a non-PD one", {
  pd_matrix <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
  result_ok <- assess_single_cholesky_decomposition(pd_matrix, list())
  expect_true(result_ok$decomposition_successful)
  expect_true(result_ok$reconstruction_error < 1e-8)

  non_pd_matrix <- matrix(c(1, 2, 2, 1), nrow = 2)
  result_fail <- assess_single_cholesky_decomposition(non_pd_matrix, list())
  expect_false(result_fail$decomposition_successful)
})

test_that("assess_depth_correlation_preservation() compares a depth-specific matrix against the original global correlation matrix", {
  original_correlations <- list(global_correlation_matrix = matrix(
    c(1, 0.5, 0.5, 1), nrow = 2, dimnames = list(c("a", "b"), c("a", "b"))
  ))
  depth_cor <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  result <- assess_depth_correlation_preservation(depth_cor, original_correlations, depth = 10, criteria = list())
  expect_equal(result$correlation_differences, 0, tolerance = 1e-8)
  expect_true(result$depth_specific_assessment)
})

test_that("assess_correlation_stability_across_depths() computes a real stability score across depth bins", {
  set.seed(1)
  n <- 100
  sim_data <- data.frame(
    hzdept_r = rep(seq(0, 90, by = 10), each = 10),
    a = stats::rnorm(n),
    b = stats::rnorm(n)
  )
  sim_data$c <- sim_data$a + stats::rnorm(n, sd = 0.1)
  result <- assess_correlation_stability_across_depths(sim_data, c("a", "b", "c"), list())
  expect_true(is.finite(result$stability_score))
})

## ----------------------------------------------------------------------
## 4c. GP performance cluster
## ----------------------------------------------------------------------

make_diag_gp_training_df <- function() {
  data.frame(
    cokey = rep(c("1", "2", "3"), each = 4),
    hzdept_r = rep(c(0, 20, 50, 100), 3),
    clay_pct = c(10, 15, 22, 30, 12, 16, 24, 32, 8, 14, 20, 28),
    stringsAsFactors = FALSE
  )
}

test_that("assess_single_gp_performance()/calculate_overall_gp_performance() derive real RMSE/R-squared from a fitted GP model", {
  fit_result <- fit_individual_gp_model(make_diag_gp_training_df(), "clay_pct", optimize_hyperparameters = FALSE)
  perf <- assess_single_gp_performance(fit_result$model, list())
  expect_true(is.finite(perf$training_rmse))
  expect_true(is.finite(perf$r_squared))

  overall <- calculate_overall_gp_performance(list(clay_pct = list(group1 = perf)))
  expect_equal(overall$mean_r_squared, perf$r_squared, tolerance = 1e-8)
})

test_that("assess_trend_realism()/summarize_constraint_violations() flag unrealistic predicted trends", {
  realistic_predictions <- c(10, 15, 20, 25, 30)
  result_ok <- assess_trend_realism(realistic_predictions, c(0, 20, 50, 100, 150), "clay_pct", list())
  expect_true(result_ok$realistic_trend)

  unrealistic_predictions <- c(10, 150, 5, 200, -50)  # noisy AND out of clay_pct's 0-100 range
  result_bad <- assess_trend_realism(unrealistic_predictions, c(0, 20, 50, 100, 150), "clay_pct", list())
  expect_false(result_bad$realistic_trend)

  summary_result <- summarize_constraint_violations(list(clay_pct = list(group1 = result_ok, group2 = result_bad)))
  expect_equal(summary_result$overall_compliance, 0.5)
  expect_true("clay_pct" %in% summary_result$violation_types)
})

test_that("validate_single_gp_predictions() returns a bounded smoothness score for a fitted GP model's predictions", {
  fit_result <- fit_individual_gp_model(make_diag_gp_training_df(), "clay_pct", optimize_hyperparameters = FALSE)
  result <- validate_single_gp_predictions(fit_result$model, list(test_depths = c(0, 20, 50, 100)))
  expect_true(is.finite(result$prediction_smoothness))
  expect_true(result$prediction_smoothness >= 0 && result$prediction_smoothness <= 1)
})

## ----------------------------------------------------------------------
## 4d. Cross-property / pedological realism cluster
## ----------------------------------------------------------------------

test_that("assess_cross_property_constraints() detects texture sum-to-100 violations and honors explicit relationship rules", {
  sim_data <- data.frame(sandtotal = c(40, 40), claytotal = c(20, 20), silttotal = c(40, 10))  # row 2 sums to 70
  result <- assess_cross_property_constraints(sim_data, cross_property_rules = list())
  expect_equal(result$texture_sum_violations, 1)

  rules <- list(list(properties = c("sandtotal", "claytotal", "silttotal"), type = "sum",
                     expected_sum = 100, tolerance = 0.05))
  result_with_rule <- assess_cross_property_constraints(sim_data, rules)
  expect_equal(result_with_rule$relationship_violations, 1)
})

test_that("detect_distribution_anomalies() counts IQR outliers across numeric columns", {
  set.seed(1)
  sim_data <- data.frame(a = c(stats::rnorm(95, 0, 1), rep(1000, 5)))
  result <- detect_distribution_anomalies(sim_data, distribution_checks = NULL)
  expect_true(result$anomalies_detected > 0)
  expect_true("outliers" %in% result$anomaly_types)
})

test_that("validate_horizon_characteristics() scores surface/subsurface horizons against realistic ranges", {
  # get_realistic_property_ranges() keys use its own naming convention
  # ("ph", "clay_total", ...), distinct from the package's usual
  # "ph1to1h2o"/"claytotal" column names elsewhere - a pre-existing naming
  # inconsistency, not something to paper over here. Use its actual key
  # names so this test exercises the range-matching logic itself.
  sim_data <- data.frame(hzdept_r = c(5, 5, 50, 50), ph = c(6.5, 6.5, 6.5, 6.5))
  result <- validate_horizon_characteristics(sim_data, list())
  expect_true(is.finite(result$surface_horizon_quality))
  expect_true(is.finite(result$subsurface_quality))
})

test_that("assess_pedological_relationships() rewards the expected clay-CEC and OM-depth directions", {
  sim_data <- data.frame(
    claytotal = c(10, 20, 30, 40, 50),
    cec7 = c(5, 10, 15, 20, 25),
    om = c(5, 4, 3, 2, 1),
    hzdept_r = c(0, 20, 40, 60, 80)
  )
  result <- assess_pedological_relationships(sim_data, list())
  expect_true(result$chemical_relationship_quality > 0.5)
  expect_true(result$physical_relationship_quality > 0.5)
})

test_that("validate_simulation_depth_trends() rewards simulated depth trends that mirror the original data", {
  set.seed(1)
  original_data <- data.frame(hzdept_r = rep(c(0, 20, 40, 60, 80), each = 4),
                              claytotal = rep(c(10, 20, 30, 40, 50), each = 4) + stats::rnorm(20, 0, 0.5))
  simulation_data <- original_data
  simulation_data$claytotal <- simulation_data$claytotal + stats::rnorm(20, 0, 0.5)

  result <- validate_simulation_depth_trends(simulation_data, original_data, list())
  expect_true(result$depth_trend_realism > 0.7)
})

## ----------------------------------------------------------------------
## 4e. Workflow scoring cluster
## ----------------------------------------------------------------------

test_that("calculate_workflow_performance() derives efficiency/coverage from real components/validation_results content", {
  components <- list(monte_carlo_results = list(x = 1), correlation_matrices = NULL, gp_models = list(y = 1))
  validation_results <- list(
    monte_carlo_validation = list(passed = TRUE),
    correlation_validation = list(validation_skipped = TRUE),
    gp_models_validation = list(validation_failed = TRUE)
  )
  result <- calculate_workflow_performance(components, validation_results)
  expect_equal(result$processing_efficiency, 2 / 3, tolerance = 1e-8)
  expect_equal(result$validation_coverage, 1 / 3, tolerance = 1e-8)
})

test_that("perform_gp_cross_validation() runs real k-fold CV per property against training_data", {
  set.seed(1)
  training_data <- data.frame(
    hzdept_r = seq(0, 100, length.out = 15),
    clay_pct = seq(10, 30, length.out = 15) + stats::rnorm(15, 0, 0.2)
  )
  gp_models <- list(clay_pct = list(type = "stratified_grouped", models = list(group1 = list())))
  result <- perform_gp_cross_validation(gp_models, training_data, list(n_folds = 3))
  expect_true(is.finite(result$cv_rmse))
  expect_true(is.finite(result$cv_r_squared))
})

## ----------------------------------------------------------------------
## 4f. Report generation cluster
## ----------------------------------------------------------------------

test_that("generate_markdown_report() writes a real, non-empty markdown file", {
  report_content <- list(executive_summary = list(overall_quality = 0.9), recommendations = c("Do X", "Do Y"))
  tmp_file <- tempfile(fileext = ".md")
  on.exit(unlink(tmp_file))
  result <- generate_markdown_report(report_content, tmp_file)
  expect_true(result)
  expect_true(file.exists(tmp_file))
  lines <- readLines(tmp_file)
  expect_true(any(grepl("executive summary", lines, ignore.case = TRUE)))
})

test_that("generate_html_report() writes a real, non-empty HTML file", {
  report_content <- list(executive_summary = list(overall_quality = 0.9))
  tmp_file <- tempfile(fileext = ".html")
  on.exit(unlink(tmp_file))
  result <- generate_html_report(report_content, tmp_file)
  expect_true(result)
  expect_true(file.exists(tmp_file))
  content <- paste(readLines(tmp_file), collapse = "\n")
  expect_true(grepl("<html", content, ignore.case = TRUE))
})

test_that("generate_pdf_report() writes a real, non-empty PDF file", {
  report_content <- list(executive_summary = list(overall_quality = 0.9))
  tmp_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp_file))
  result <- generate_pdf_report(report_content, tmp_file)
  expect_true(result)
  expect_true(file.exists(tmp_file))
  expect_true(file.info(tmp_file)$size > 0)
})

test_that("render_value_as_lines() renders nested lists without calling print() on non-atomic objects", {
  weird_object <- structure(list(), class = "some_weird_class")
  lines <- render_value_as_lines(list(a = 1, b = list(c = weird_object)))
  expect_true(any(grepl("some_weird_class", lines)))
})

test_that("theme_soil_diagnostics()/save_diagnostic_plot() produce a real usable ggplot theme and handle NULL output_dir", {
  skip_if_not_installed("ggplot2")
  theme_obj <- theme_soil_diagnostics()
  expect_s3_class(theme_obj, "theme")

  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_identical(save_diagnostic_plot(p, "unused.png", NULL), p)
})

test_that("generate_comprehensive_diagnostics() degrades each plot category to an empty list when its inputs are missing", {
  result <- generate_comprehensive_diagnostics(components = list(), validation_results = list(), output_dir = NULL)
  expect_length(result$monte_carlo_plots, 0)
  expect_length(result$gp_model_plots, 0)
  expect_length(result$soil_science_plots, 0)
  expect_length(result$summary_plots, 0)
  expect_equal(result$plot_generation_status, "complete")
})

test_that("generate_comprehensive_diagnostics() produces a real Monte Carlo distribution plot when simulation_data is available", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  components <- list(simulation_data = data.frame(
    cokey = "1", simulation_number = 1:20,
    clay_pct = stats::rnorm(20, 20, 3), sand_pct = stats::rnorm(20, 40, 3)
  ))
  result <- generate_comprehensive_diagnostics(components, validation_results = list(), output_dir = NULL)
  expect_true(!is.null(result$monte_carlo_plots$property_distributions))
  expect_s3_class(result$monte_carlo_plots$property_distributions, "ggplot")
})

test_that("create_quality_assessment() combines score/grade/critical-issues into the documented shape", {
  assessment <- create_quality_assessment(
    overall_score = 0.85, quality_grade = "Good",
    quality_scores = list(monte_carlo = 0.9), critical_issues = character(0),
    validation_config = list(minimum_acceptable_score = 0.7)
  )
  expect_equal(assessment$workflow_status, "PASSED")
  expect_equal(assessment$quality_grade, "Good")

  failing <- create_quality_assessment(
    overall_score = 0.5, quality_grade = "Poor",
    quality_scores = list(), critical_issues = "some issue",
    validation_config = list(minimum_acceptable_score = 0.7)
  )
  expect_equal(failing$workflow_status, "FAILED")
})
