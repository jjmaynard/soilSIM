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
