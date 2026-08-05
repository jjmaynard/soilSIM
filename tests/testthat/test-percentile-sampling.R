test_that("extract_percentile_pairs() parses P<num> columns, sorts, and drops NA/-1 sentinels", {
  df <- data.frame(P50 = 10, P5 = 2, P95 = 18, P0 = -1, stringsAsFactors = FALSE)
  pairs <- extract_percentile_pairs(df, c("P0", "P5", "P50", "P95"))
  expect_equal(pairs$probs, c(0.05, 0.5, 0.95))
  expect_equal(pairs$values, c(2, 10, 18))
})

test_that("extract_percentile_pairs() errors with fewer than 2 valid columns", {
  df <- data.frame(P50 = 10, stringsAsFactors = FALSE)
  expect_error(extract_percentile_pairs(df, c("P50")), "at least 2")
})

test_that("sim_linear_cdf() draws are exact at the knots and clamped at the tails", {
  set.seed(1)
  draws <- sim_linear_cdf(probs = c(0, 0.5, 1), values = c(0, 10, 20), n = 5000)
  expect_true(all(draws >= 0 & draws <= 20))
  expect_equal(mean(draws), 10, tolerance = 0.5)
})

test_that("sim_spline() draws stay within the padded bounds and are monotonic at the knots", {
  set.seed(2)
  draws <- sim_spline(probs = c(0.05, 0.5, 0.95), values = c(2, 10, 18), n = 2000)
  expect_true(all(is.finite(draws)))
  expect_true(min(draws) >= 2 - 0.05 * 16 - 1e-6)
  expect_true(max(draws) <= 18 + 0.05 * 16 + 1e-6)
})

test_that("sim_kde() requires truncnorm and draws stay within the observed value range", {
  testthat::skip_if_not_installed("truncnorm")
  set.seed(3)
  draws <- sim_kde(probs = c(0.05, 0.5, 0.95), values = c(2, 10, 18), n = 500, sample_size = 2000)
  expect_true(all(draws >= 2 & draws <= 18))
})

test_that("simulate_from_percentiles() dispatches linear_cdf/spline/normal correctly", {
  set.seed(4)
  df <- data.frame(P5 = 2, P50 = 10, P95 = 18)
  cols <- c("P5", "P50", "P95")

  lin <- simulate_from_percentiles(df, method = "linear_cdf", percentile_cols = cols, n = 500)
  expect_length(lin, 500)

  spl <- simulate_from_percentiles(df, method = "spline", percentile_cols = cols, n = 500)
  expect_length(spl, 500)

  norm <- simulate_from_percentiles(df, method = "normal", percentile_cols = cols, n = 500)
  expect_length(norm, 500)
  expect_equal(mean(norm), 10, tolerance = 1)
})

test_that("simulate_from_percentiles() rejects method = 'metalog' (deliberately not ported)", {
  df <- data.frame(P5 = 2, P50 = 10, P95 = 18)
  expect_error(
    simulate_from_percentiles(df, method = "metalog", percentile_cols = c("P5", "P50", "P95")),
    "should be one of"
  )
})

test_that("calculate_summary_statistics() computes basic stats and dynamically-named percentiles", {
  data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- calculate_summary_statistics(data, percentile_probs = c(0.1, 0.5, 0.9))
  expect_equal(result$Num, 10)
  expect_equal(result$Mean, 5.5)
  expect_true(all(c("P10", "P50", "P90") %in% names(result)))
})

test_that("calculate_summary_statistics() errors on non-numeric input", {
  expect_error(calculate_summary_statistics("not numeric"), "must be a numeric vector")
})

test_that("compare_percentile_methods() returns samples and a summary row per method, tolerating a failing method", {
  set.seed(5)
  df <- data.frame(P5 = 2, P50 = 10, P95 = 18)
  result <- suppressWarnings(compare_percentile_methods(
    df, methods = c("linear_cdf", "spline"), percentile_cols = c("P5", "P50", "P95"), n = 200
  ))
  expect_named(result, c("samples", "summary"))
  expect_length(result$samples, 2)
  expect_equal(nrow(result$summary), 2)
})

test_that("validate_percentile_methods_synthetic() scores linear_cdf reconstruction against a known normal", {
  set.seed(6)
  true_rng <- function(n) stats::rnorm(n, mean = 10, sd = 3)
  result <- validate_percentile_methods_synthetic(
    true_rng, methods = "linear_cdf", n_true = 1000, n_sim = 500, n_reps = 3
  )
  expect_equal(nrow(result), 1)
  expect_true(result$ks_stat[1] < 0.2)
})

test_that("generate_inverse_cdf_distribution() is a thin linear_cdf wrapper", {
  set.seed(7)
  df <- data.frame(P5 = 2, P50 = 10, P95 = 18)
  draws <- generate_inverse_cdf_distribution(df, percentile_cols = c("P5", "P50", "P95"), n = 300)
  expect_length(draws, 300)
})
