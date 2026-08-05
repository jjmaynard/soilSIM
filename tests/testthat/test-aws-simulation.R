test_that("van_genuchten() returns theta_s at h = 0 regardless of alpha/n", {
  expect_equal(van_genuchten(h = 0, alpha = 0.02, n = 1.3, theta_r = 0.05, theta_s = 0.45), 0.45)
  expect_equal(van_genuchten(h = 0, alpha = 0.5, n = 2.5, theta_r = 0.1, theta_s = 0.4), 0.4)
})

test_that("van_genuchten() approaches theta_r as h becomes very negative", {
  theta_r <- 0.05
  theta_s <- 0.45
  wet <- van_genuchten(h = -1, alpha = 0.02, n = 1.3, theta_r = theta_r, theta_s = theta_s)
  dry <- van_genuchten(h = -1e15, alpha = 0.02, n = 1.3, theta_r = theta_r, theta_s = theta_s)
  expect_true(dry < wet)
  expect_equal(dry, theta_r, tolerance = 1e-3)
})

test_that("van_genuchten() is monotonically non-increasing in |h|", {
  h_seq <- -c(1, 10, 100, 1000, 10000)
  theta <- van_genuchten(h = h_seq, alpha = 0.02, n = 1.3, theta_r = 0.05, theta_s = 0.45)
  expect_true(all(diff(theta) <= 0))
})

test_that("simulate_vg_aws() returns a named list with AWHC = theta_fc - theta_pwp for each row", {
  data <- make_rosetta_shaped_data()
  result <- simulate_vg_aws(data, n_simulations = 20)

  expect_type(result, "list")
  expect_length(result, nrow(data))
  expect_true(all(grepl("^testseries_0_|^testseries_20_", names(result))))

  for (df in result) {
    expect_equal(nrow(df), 20)
    expect_equal(df$AWHC, df$theta_fc - df$theta_pwp)
    expect_true(all(c("alpha", "n", "theta_r", "theta_s", "sim_num", "theta_fc", "theta_pwp", "AWHC") %in% names(df)))
  }
})

test_that("simulate_vg_aws() skips rows with any missing van Genuchten parameter", {
  data <- make_rosetta_shaped_data()
  data$theta_s[1] <- NA
  result <- simulate_vg_aws(data, n_simulations = 5)
  expect_length(result, 1)
  expect_true(grepl("^testseries_20_", names(result)))
})

test_that("calculate_aws_df() requires the live ROSETTA (handbook60.org) service", {
  testthat::skip_if_offline()
  sim_data_df <- make_aws_texture_data()
  result <- calculate_aws_df(sim_data_df)

  # Long format: one row per cokey per depth slab actually spanned by the
  # component's horizons (not one row per cokey - see calculate_aws_df()'s
  # @return doc), with a single AWHC value column.
  expect_s3_class(result, "data.frame")
  expect_setequal(names(result), c("cokey", "top", "bottom", "AWHC"))
  expect_true(nrow(result) >= 1)
  expect_true(all(result$cokey == "1"))
  expect_true(all(result$AWHC >= 0 & result$AWHC <= 1))
})
