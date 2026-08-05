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

test_that("simulate_cokey_generalized() skips a row with no recognized properties rather than erroring", {
  sim_cokey <- data.frame(
    genhz = "A", sim_comppct = 10, compname = "x", mukey = "1", cokey = "1",
    hzdept_r = 0, hzdepb_r = 20, stringsAsFactors = FALSE
  )
  corr <- make_property_correlation_matrices()
  expect_message(result <- simulate_cokey_generalized(sim_cokey, corr), "No recognized properties")
  expect_null(result)
})
