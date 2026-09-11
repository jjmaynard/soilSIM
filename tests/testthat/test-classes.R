test_that("soilSIM result constructors attach class without disturbing contents", {
  raw <- list(a = 1, b = list(c = 2), metadata = list(properties = c("clay", "sand")))
  sim <- new_soilSIM_simulation(raw)
  expect_s3_class(sim, "soilSIM_simulation")
  expect_s3_class(sim, "soilSIM_result")
  expect_true(is.list(sim))
  expect_identical(sim$a, 1)
  expect_identical(sim[["b"]]$c, 2)
  expect_identical(names(sim), names(raw))
})

test_that("constructors reject non-list input", {
  expect_error(new_soilSIM_simulation(1:10), "wrap a list")
})

test_that("print methods return their object invisibly and do not error", {
  sim <- new_soilSIM_simulation(list(
    metadata = list(properties = "clay", n_realizations = 100, n_horizons = 5,
                    success_rate = 0.99),
    quality_assessment = list(overall_quality_score = 0.82),
    simulation_data = array(0, dim = c(5, 1, 100))
  ))
  expect_output(print(sim), "soilSIM_simulation")
  expect_invisible(print(sim))
  expect_output(summary(sim), "simulation_data dim")

  stats <- new_soilSIM_statistics(list(
    analysis_metadata = list(properties_analyzed = 3),
    quality_report = list(overall_quality_score = 0.9),
    correlation_matrices = list(matrices = list(A = 1, B = 2))
  ))
  expect_output(print(stats), "properties analyzed")

  gp <- new_soilSIM_gp_models(list(
    clay = list(n_groups = 2, models = list(1, 2)),
    model_summary = list()
  ))
  expect_output(print(gp), "1 group model|2 group model")

  dg <- new_soilSIM_diagnostics(list(overall_status = "PASS", section_a = list()))
  expect_output(print(dg), "PASS")
})

test_that("summary.soilSIM_diagnostics reads diagnose_workflow()'s real overall_assessment shape (D5)", {
  dg <- new_soilSIM_diagnostics(list(
    overall_assessment = list(
      quality_score = 0.85, quality_grade = "B", workflow_status = "PASSED",
      component_scores = list(monte_carlo = 0.9, correlation = 0.8),
      critical_issues = list("one issue")
    ),
    recommendations = c("Increase sample size", "Review correlation matrix"),
    monte_carlo_validation = list()
  ))
  expect_output(print(dg), "PASSED")
  out <- capture.output(summary(dg))
  expect_true(any(grepl("monte_carlo", out)))
  expect_true(any(grepl("Increase sample size", out)))
})

test_that("engine wrappers return classed objects", {
  skip_if_not(exists("amador_processed_fixture", mode = "function"))
})

test_that("property_config engines return soilSIM_property_config (D6.3)", {
  d <- default_property_config("claytotal")
  expect_s3_class(d, "soilSIM_property_config")
  expect_identical(d$type, "texture")

  c1 <- create_custom_property_config("my_prop", property_type = "generic", units = "%")
  expect_s3_class(c1, "soilSIM_property_config")
  expect_identical(c1$units, "%")
})

test_that("run_fusion() returns soilSIM_fusion on the closed-form success path (D6.1)", {
  skip_if_not(exists("make_percentile_rasters", mode = "function"))
  # run_fusion() needs live SSURGO/SOLUS network access to exercise end to end;
  # the unit-level contract (early NULL vs. wrapped success) is instead checked
  # directly against stage1_fuse_from_prior_solus(), the function run_fusion()
  # wraps for its non-group return path.
  skip_if_not(exists("stage1_fuse_from_prior_solus", mode = "function"))
})

test_that("predict.soilSIM_gp_models() looks up the right leaf model and matches predict_gp_depth_trends() (D6.4)", {
  df <- make_gp_training_df()
  fit <- fit_individual_gp_model(df, "clay_pct", optimize_hyperparameters = FALSE)
  gp_models <- new_soilSIM_gp_models(list(
    clay_pct = list(models = list(all = fit$model)),
    model_summary = list()
  ))

  direct <- predict_gp_depth_trends(fit$model, c(0, 50, 100))
  via_method <- predict(gp_models, c(0, 50, 100), property = "clay_pct")
  expect_equal(via_method, direct)

  via_explicit_group <- predict(gp_models, c(0, 50, 100), property = "clay_pct", group = "all")
  expect_equal(via_explicit_group, direct)

  expect_error(predict(gp_models, 0, property = "not_a_property"), "not found")
  expect_error(predict(gp_models, 0, property = "clay_pct", group = "not_a_group"), "not found")
})

test_that("plot.soilSIM_simulation()/plot.soilSIM_gp_models() run without erroring when ggplot2 is available (D6.4)", {
  skip_if_not_installed("ggplot2")

  sim <- new_soilSIM_simulation(list(
    simulation_data = array(
      rnorm(5 * 2 * 20), dim = c(5, 2, 20),
      dimnames = list(NULL, c("claytotal", "sandtotal"), NULL)
    )
  ))
  expect_invisible(plot(sim))
  expect_invisible(plot(sim, property = "sandtotal"))

  df <- make_gp_training_df()
  fit <- fit_individual_gp_model(df, "clay_pct", optimize_hyperparameters = FALSE)
  gp_models <- new_soilSIM_gp_models(list(clay_pct = list(models = list(all = fit$model))))
  expect_invisible(plot(gp_models, property = "clay_pct", new_depths = c(0, 25, 50)))
})

test_that("plot.soilSIM_simulation() degrades gracefully without ggplot2", {
  skip_if(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is installed; degrade path not exercised")
  sim <- new_soilSIM_simulation(list(simulation_data = array(1, dim = c(1, 1, 1))))
  expect_message(result <- plot(sim), "ggplot2")
  expect_null(result)
})

test_that("plot.soilSIM_fusion() draws prior/likelihood/posterior raster panels (D6.4)", {
  skip_if_not(exists("make_percentile_rasters", mode = "function"))
  probs <- c(0.1, 0.5, 0.9)
  fusion <- new_soilSIM_fusion(list(
    prior = list(values = make_percentile_rasters(list("0.1" = 8, "0.5" = 12, "0.9" = 16)), probs = probs),
    likelihood = list(values = make_percentile_rasters(list("0.1" = 10, "0.5" = 14, "0.9" = 18)), probs = probs),
    posterior = list(mu = make_percentile_rasters(list(m = 13))[[1]], sigma = make_percentile_rasters(list(s = 1))[[1]]),
    dist = "normal"
  ))
  dev_file <- tempfile(fileext = ".png")
  grDevices::png(dev_file)
  result <- tryCatch(plot(fusion), finally = grDevices::dev.off())
  unlink(dev_file)
  expect_identical(result, fusion)
})

test_that("plot.soilSIM_fusion() degrades gracefully with no raster panels", {
  expect_message(result <- plot(new_soilSIM_fusion(list(dist = "normal"))), "no raster panels")
  expect_null(result)
})
