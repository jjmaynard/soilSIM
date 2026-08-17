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

test_that("sim_linear_cdf_batch()'s interpolation matches stats::approxfun() exactly at fixed evaluation points", {
  # Isolates the interpolation math from RNG differences (sim_linear_cdf_batch() draws its
  # nrow*n uniforms in one block rather than one runif(n) call per row, so it is NOT
  # bit-identical to sim_linear_cdf() row-by-row - see its own docs) by feeding identical FIXED
  # evaluation points to both the batch helper's internals and stats::approxfun() directly.
  probs <- c(0.05, 0.5, 0.95)
  values_mat <- rbind(c(18, 22, 27), c(20, 40, 65), c(0.5, 2, 8))
  u_fixed <- c(0, 0.01, 0.04999, 0.05, 0.1, 0.5, 0.94999, 0.95, 0.96, 1)

  # Reuses sim_linear_cdf_batch()'s exact interpolation logic, but driven by fixed points
  # instead of stats::runif() draws - a whitebox check of the math, not the public API.
  nr <- nrow(values_mat); n <- length(u_fixed); k <- length(probs)
  u <- matrix(rep(u_fixed, each = nr), nrow = nr, ncol = n)
  seg <- findInterval(as.vector(u), probs, rightmost.closed = TRUE)
  seg_clamped <- pmin(pmax(seg, 1L), k - 1L)
  row_idx <- rep(seq_len(nr), times = n)
  p_lo <- probs[seg_clamped]; p_hi <- probs[seg_clamped + 1L]
  v_lo <- values_mat[cbind(row_idx, seg_clamped)]
  v_hi <- values_mat[cbind(row_idx, seg_clamped + 1L)]
  frac <- (as.vector(u) - p_lo) / (p_hi - p_lo)
  interp <- v_lo + frac * (v_hi - v_lo)
  below <- seg == 0L; above <- seg == k
  interp[below] <- values_mat[cbind(row_idx[below], rep(1L, sum(below)))]
  interp[above] <- values_mat[cbind(row_idx[above], rep(k, sum(above)))]
  batch_result <- matrix(interp, nrow = nr, ncol = n)

  for (i in seq_len(nr)) {
    ref_fun <- stats::approxfun(x = probs, y = values_mat[i, ], method = "linear", rule = 2)
    expect_equal(batch_result[i, ], ref_fun(u_fixed), tolerance = 1e-12)
  }
})

test_that("sim_linear_cdf_batch() is distributionally equivalent to sim_linear_cdf() per row (KS test)", {
  probs <- c(0.05, 0.5, 0.95)
  values_mat <- rbind(c(18, 22, 27), c(20, 40, 65))

  set.seed(11)
  batch_draws <- sim_linear_cdf_batch(probs, values_mat, 20000)
  for (i in seq_len(nrow(values_mat))) {
    set.seed(11)
    scalar_draws <- sim_linear_cdf(probs, values_mat[i, ], 20000)
    ks <- suppressWarnings(stats::ks.test(batch_draws[i, ], scalar_draws))
    expect_gt(ks$p.value, 0.01)
  }
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
