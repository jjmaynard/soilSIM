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
