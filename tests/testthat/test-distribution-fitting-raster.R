test_that("fit_normal_raster()/quantile_normal_raster() match the scalar fit_normal_triplet()/quantile_normal() cell-by-cell", {
  rasters <- make_percentile_rasters(c(P5 = 2, P50 = 10, P95 = 18))
  fit_r <- fit_normal_raster(rasters$P5, rasters$P50, rasters$P95, 0.05, 0.95)
  fit_scalar <- fit_normal_triplet(2, 10, 18, 0.05, 0.95)

  expect_equal(unique(terra::values(fit_r$mu))[1], fit_scalar$mean)
  expect_equal(unique(terra::values(fit_r$sigma))[1], fit_scalar$sd)

  q_r <- quantile_normal_raster(fit_r, 0.75)
  q_scalar <- quantile_normal(fit_scalar, 0.75)
  expect_equal(unique(terra::values(q_r))[1], unname(q_scalar))
})

test_that("quantile_linear_cdf_raster() matches the scalar quantile_linear_cdf() cell-by-cell", {
  rasters <- make_percentile_rasters(c(P0 = 0, P50 = 10, P100 = 20))
  q_r <- quantile_linear_cdf_raster(list(rasters$P0, rasters$P50, rasters$P100), c(0, 0.5, 1), 0.25)
  q_scalar <- quantile_linear_cdf(c(0, 0.5, 1), c(0, 10, 20), 0.25)
  expect_equal(unique(terra::values(q_r))[1], q_scalar)
})

test_that("fit_beta_mom_raster() matches the scalar fit_beta_mom() cell-by-cell", {
  mean_r <- make_percentile_rasters(c(m = 0.4))$m
  var_r <- make_percentile_rasters(c(v = 0.02))$v
  fit_r <- fit_beta_mom_raster(mean_r, var_r)
  fit_scalar <- fit_beta_mom(0.4, 0.02)
  expect_equal(unique(terra::values(fit_r$alpha))[1], fit_scalar$shape1)
  expect_equal(unique(terra::values(fit_r$beta))[1], fit_scalar$shape2)
})

test_that("raster_digamma()/raster_trigamma() match base digamma()/trigamma()", {
  r <- make_percentile_rasters(c(x = 3.5))$x
  expect_equal(unique(terra::values(raster_digamma(r)))[1], digamma(3.5))
  expect_equal(unique(terra::values(raster_trigamma(r)))[1], trigamma(3.5))
})

test_that("fit_beta_mle_newton_raster()/quantile_beta_mle_newton_raster() match the scalar fit_beta_mle_newton()/quantile_beta() cell-by-cell", {
  rasters <- make_percentile_rasters(c(P5 = 20, P50 = 40, P95 = 65))
  bounds <- c(0, 100)
  fit_r <- fit_beta_mle_newton_raster(list(rasters$P5, rasters$P50, rasters$P95), bounds)
  fit_scalar <- fit_beta_mle_newton(c(20, 40, 65), bounds)

  expect_equal(unique(terra::values(fit_r$alpha))[1], fit_scalar$shape1, tolerance = 1e-6)
  expect_equal(unique(terra::values(fit_r$beta))[1], fit_scalar$shape2, tolerance = 1e-6)

  q_r <- quantile_beta_mle_newton_raster(fit_r, 0.5)
  q_scalar <- qbeta(0.5, fit_scalar$shape1, fit_scalar$shape2)
  expect_equal(unique(terra::values(q_r))[1], q_scalar, tolerance = 1e-6)
})

test_that("fit_metalog_linear_raster()/quantile_metalog_linear_raster() match the scalar fit_metalog_linear()/quantile_metalog_linear() cell-by-cell", {
  rasters <- make_percentile_rasters(c(P10 = 3, P50 = 10, P90 = 22))
  probs <- c(0.1, 0.5, 0.9)
  bounds <- c(0, 100)
  fit_r <- fit_metalog_linear_raster(list(rasters$P10, rasters$P50, rasters$P90), probs, bounds, "b")
  fit_scalar <- fit_metalog_linear(c(3, 10, 22), probs, bounds, "b")

  q_r <- quantile_metalog_linear_raster(fit_r, 0.5, bounds, "b")
  q_scalar <- quantile_metalog_linear(fit_scalar, 0.5)
  expect_equal(unique(terra::values(q_r))[1], unname(q_scalar), tolerance = 1e-6)
})

test_that("check_metalog_feasibility_raster() flags no infeasible cells for a well-behaved fit", {
  rasters <- make_percentile_rasters(c(P10 = 3, P50 = 10, P90 = 22))
  probs <- c(0.1, 0.5, 0.9)
  bounds <- c(0, 100)
  fit_r <- fit_metalog_linear_raster(list(rasters$P10, rasters$P50, rasters$P90), probs, bounds, "b")
  infeasible_r <- check_metalog_feasibility_raster(fit_r, bounds, "b")
  expect_false(any(terra::values(infeasible_r)))
})

test_that("quantile_metalog_linear_with_fallback() uses the metalog fit when feasible", {
  rasters <- make_percentile_rasters(c(P10 = 3, P50 = 10, P90 = 22))
  probs <- c(0.1, 0.5, 0.9)
  bounds <- c(0, 100)
  fit_r <- fit_metalog_linear_raster(list(rasters$P10, rasters$P50, rasters$P90), probs, bounds, "b")
  infeasible_r <- check_metalog_feasibility_raster(fit_r, bounds, "b")

  result <- quantile_metalog_linear_with_fallback(
    fit_r, infeasible_r, list(rasters$P10, rasters$P50, rasters$P90), probs, 0.5, bounds, "b"
  )
  direct <- quantile_metalog_linear_raster(fit_r, 0.5, bounds, "b")
  expect_equal(terra::values(result$value), terra::values(direct))
  expect_false(any(terra::values(result$used_fallback)))
})
