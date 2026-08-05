test_that("fit_normal_triplet/quantile_normal exact closed form matches hand computation", {
  # Symmetric-around-the-median triplet (5, 10, 15) is internally consistent
  # with an actual Normal distribution, so the fit round-trips exactly at all
  # three input percentiles - an asymmetric triplet (e.g. 5/10/20) is NOT
  # generally consistent with any real Normal and would not round-trip at
  # the tails (the closed-form fit is a best univariate summary, not an
  # interpolant, for inputs that aren't actually Normal-consistent).
  fit <- fit_normal_triplet(5, 10, 15, 0.05, 0.95)
  expect_equal(fit$mean, 10)
  expect_equal(fit$sd, (15 - 5) / (qnorm(0.95) - qnorm(0.05)))
  expect_equal(quantile_normal(fit, 0.5), 10, tolerance = 1e-9)
  expect_equal(quantile_normal(fit, 0.95), 15, tolerance = 1e-9)
  expect_equal(quantile_normal(fit, 0.05), 5, tolerance = 1e-9)
})

test_that("fit_beta_mle_newton matches fitdistrplus::fitdist on the same percentile values", {
  vals <- c(12, 28, 55)
  ours <- fit_beta_mle_newton(vals, bounds = c(0, 100))
  theirs <- fitdistrplus::fitdist(vals / 100, "beta")$estimate
  expect_equal(ours$shape1, unname(theirs[["shape1"]]), tolerance = 1e-2)
  expect_equal(ours$shape2, unname(theirs[["shape2"]]), tolerance = 1e-2)
})

test_that("fit_beta_mle_newton_vec matches fit_beta_mle_newton row-by-row", {
  mat <- rbind(c(12, 28, 55), c(5, 20, 45), c(30, 50, 70))
  vec_fit <- fit_beta_mle_newton_vec(mat, bounds = c(0, 100))
  for (i in seq_len(nrow(mat))) {
    single <- fit_beta_mle_newton(mat[i, ], bounds = c(0, 100))
    expect_equal(vec_fit$shape1[i], single$shape1, tolerance = 1e-9)
    expect_equal(vec_fit$shape2[i], single$shape2, tolerance = 1e-9)
  }
})

test_that("metalog fit round-trips exactly at its own fitting probabilities", {
  fit <- fit_metalog_linear(c(5, 20, 45), c(0.05, 0.5, 0.95), bounds = c(0, 100), boundedness = "b")
  recovered <- quantile_metalog_linear(fit, c(0.05, 0.5, 0.95))
  expect_equal(recovered, c(5, 20, 45), tolerance = 1e-6)
})

test_that("check_metalog_feasible flags a constructed non-monotonic case, and fallback matches linear_cdf exactly", {
  fit <- fit_metalog_linear(c(0.1, 0.5, 99.9), c(0.05, 0.5, 0.95), bounds = c(0, 100), boundedness = "b")
  expect_true(check_metalog_feasible(fit))

  q <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  fb <- quantile_metalog_with_fallback(fit, infeasible = TRUE, full_probs = c(0.05, 0.5, 0.95), full_values = c(0.1, 0.5, 99.9), q = q)
  lc <- quantile_linear_cdf(c(0.05, 0.5, 0.95), c(0.1, 0.5, 99.9), q)
  expect_equal(fb, lc)
})

test_that("check_metalog_feasible does not flag a well-behaved case", {
  fit <- fit_metalog_linear(c(5, 20, 45), c(0.05, 0.5, 0.95), bounds = c(0, 100), boundedness = "b")
  expect_false(check_metalog_feasible(fit))
})

test_that("fit_percentile_triplet dispatches every family and returns valid=FALSE on non-finite input", {
  fam_normal <- fit_percentile_triplet(5, 10, 20, "normal")
  expect_equal(fam_normal$family, "normal")
  expect_true(fam_normal$valid)

  fam_beta <- fit_percentile_triplet(12, 28, 55, "beta", bounds = c(0, 100))
  expect_equal(fam_beta$family, "beta")
  expect_true(fam_beta$fit$shape1 > 0 && fam_beta$fit$shape2 > 0)

  fam_bad <- fit_percentile_triplet(NA, 10, 20, "normal")
  expect_false(fam_bad$valid)

  fam_lognormal_guard <- fit_percentile_triplet(-1, 0, 5, "lognormal")
  expect_equal(fam_lognormal_guard$family, "normal")  # non-positive values fall back to normal
})

test_that("fit_percentile_triplet resolves 'auto' via skew proxy and never resolves to metalog", {
  symmetric <- fit_percentile_triplet(5, 10, 15, "auto")
  expect_equal(symmetric$family, "normal")

  skewed <- fit_percentile_triplet(1, 2, 100, "auto")
  expect_equal(skewed$family, "lognormal")

  bounded <- fit_percentile_triplet(12, 28, 55, "auto", bounds = c(0, 100))
  expect_equal(bounded$family, "beta")
})

test_that("fit_percentile_triplet handles a zero-width interval as a degenerate constant", {
  deg <- fit_percentile_triplet(10, 10, 10, "beta", bounds = c(0, 100))
  expect_equal(deg$family, "triangular")
  expect_equal(deg$fit, list(min = 10, mode = 10, max = 10))
  expect_equal(quantile_from_fit(c(0.1, 0.5, 0.9), deg$family, deg$fit), c(10, 10, 10))
})

test_that("quantile_from_fit returns NA (not an error) for a NULL/invalid fit", {
  expect_true(all(is.na(quantile_from_fit(c(0.1, 0.5, 0.9), "normal", NULL))))
})

test_that("validate_fit_parameters rejects invalid parameter sets and accepts valid ones", {
  expect_true(validate_fit_parameters("beta", list(shape1 = 2, shape2 = 3, lower = 0, upper = 100))$valid)
  expect_false(validate_fit_parameters("beta", list(shape1 = -1, shape2 = 3, lower = 0, upper = 100))$valid)
  expect_false(validate_fit_parameters("beta", list(shape1 = 2, shape2 = 3, lower = 100, upper = 0))$valid)
  expect_true(validate_fit_parameters("normal", list(mean = 5, sd = 2))$valid)
  expect_false(validate_fit_parameters("normal", list(mean = 5, sd = -1))$valid)
  expect_false(validate_fit_parameters("triangular", list(min = 10, mode = 20, max = 5))$valid)

  # An infeasible metalog fit is still "valid" (transparently handled via the
  # linear_cdf fallback), it just carries a diagnostic message.
  infeasible_fit <- fit_metalog_linear(c(0.1, 0.5, 99.9), c(0.05, 0.5, 0.95), bounds = c(0, 100), boundedness = "b")
  infeasible_fit$infeasible <- TRUE
  infeasible_fit$fallback_probs <- c(0.05, 0.5, 0.95)
  infeasible_fit$fallback_values <- c(0.1, 0.5, 99.9)
  result <- validate_fit_parameters("metalog", infeasible_fit)
  expect_true(result$valid)
  expect_match(result$message, "fallback")
})

test_that("ilr_forward/ilr_inverse round-trip to ~1e-12 and preserve sum/bounds", {
  set.seed(7)
  n <- 200
  clay <- runif(n, 5, 40)
  sand <- runif(n, 5, 40)
  silt <- 100 - clay - sand  # guaranteed > 0 given the ranges above, and sums to exactly 100

  z <- ilr_forward(clay, sand, silt)
  back <- ilr_inverse(z[, "z1"], z[, "z2"], total = 100)

  expect_equal(unname(back[, "clay"]), clay, tolerance = 1e-9)
  expect_equal(unname(back[, "sand"]), sand, tolerance = 1e-9)
  expect_equal(unname(back[, "silt"]), silt, tolerance = 1e-9)
  expect_equal(rowSums(back), rep(100, n), tolerance = 1e-9)
  expect_true(all(back >= 0 & back <= 100))
})

test_that("estimate_ilr_moments_mc recovers the right mean under symmetric low/high and returns a valid PD Sigma", {
  set.seed(11)
  z <- qnorm(0.95)
  moments <- estimate_ilr_moments_mc(
    low_clay = 15, rep_clay = 20, high_clay = 25,
    low_sand = 35, rep_sand = 40, high_sand = 45,
    low_silt = 35, rep_silt = 40, high_silt = 45,
    z = z, n_mc = 20000
  )
  true_z <- ilr_forward(20, 40, 40)
  expect_equal(unname(moments$mu), as.numeric(true_z), tolerance = 0.05)
  expect_true(all(eigen(moments$Sigma, only.values = TRUE)$values > 0))
  expect_equal(moments$Sigma[1, 2], moments$Sigma[2, 1])
})

test_that("sample_ilr_posterior draws sum to total and stay in bounds for every draw", {
  set.seed(13)
  mu <- ilr_forward(20, 40, 40)
  Sigma <- diag(c(0.05, 0.05))
  draws <- sample_ilr_posterior(as.numeric(mu), Sigma, n = 500, total = 100)
  expect_equal(rowSums(draws), rep(100, 500), tolerance = 1e-8)
  expect_true(all(draws >= 0 & draws <= 100))
})

test_that("ensure_positive_definite_matrix repairs a matrix with a negative eigenvalue", {
  bad <- matrix(c(1, 0.99, 0.99, 0.99, 1, 0.99, 0.99, 0.99, 1), nrow = 3)
  bad[1, 3] <- bad[3, 1] <- -0.99  # forces a negative eigenvalue while keeping unit diagonal
  repaired <- ensure_positive_definite_matrix(bad)
  expect_true(all(eigen(repaired, only.values = TRUE)$values > 0))
  expect_equal(diag(repaired), rep(1, 3), tolerance = 1e-8)
})

test_that("ensure_positive_definite_matrix preserves dimnames (previously silently dropped)", {
  m <- diag(3)
  rownames(m) <- colnames(m) <- c("a", "b", "c")
  repaired <- ensure_positive_definite_matrix(m)
  expect_equal(dimnames(repaired), list(c("a", "b", "c"), c("a", "b", "c")))
  expect_equal(repaired["a", "b"], 0, tolerance = 1e-8)
})

test_that("validate_correlation_matrix catches structural problems", {
  good <- diag(3)
  expect_true(validate_correlation_matrix(good, c("a", "b", "c"))$valid)
  expect_false(validate_correlation_matrix(matrix(1:6, nrow = 2), c("a", "b"))$valid)
  expect_false(validate_correlation_matrix(diag(3), c("a", "b"))$valid)

  non_unit_diag <- diag(3); non_unit_diag[1, 1] <- 0.5
  expect_false(validate_correlation_matrix(non_unit_diag, c("a", "b", "c"))$valid)

  non_pd <- matrix(c(1, 0.9, 0.9, 0.9, 1, 0.9, 0.9, 0.9, 1), nrow = 3)
  non_pd[1, 3] <- non_pd[3, 1] <- -0.99
  expect_false(validate_correlation_matrix(non_pd, c("a", "b", "c"))$valid)
})

test_that("estimate_correlation_matrix_robust uses empirical correlation when n is sufficient, else falls back", {
  set.seed(17)
  n <- 100
  x <- rnorm(n)
  y <- 0.8 * x + rnorm(n, sd = 0.3)
  z <- rnorm(n)
  df <- data.frame(x = x, y = y, z = z)

  result <- estimate_correlation_matrix_robust(df, min_group_n = 5)
  expect_equal(result$method, "empirical_pooled")
  expect_true(abs(result$matrix["x", "y"]) > 0.5)
  expect_true(all(eigen(result$matrix, only.values = TRUE)$values > 0))

  tiny <- df[1:3, ]
  result_small <- estimate_correlation_matrix_robust(tiny, min_group_n = 5)
  expect_equal(result_small$method, "global_fallback")
  expect_equal(result_small$matrix, diag(3), ignore_attr = TRUE)
})

test_that("estimate_correlation_matrix_robust with group_var combines qualifying groups and falls back for n=1 groups", {
  set.seed(19)
  n_per_group <- 30
  make_group <- function(rho) {
    x <- rnorm(n_per_group)
    y <- rho * x + rnorm(n_per_group, sd = sqrt(1 - rho^2))
    data.frame(x = x, y = y)
  }
  df <- rbind(
    cbind(make_group(0.7), genhz = "A"),
    cbind(make_group(0.7), genhz = "B"),
    cbind(data.frame(x = 1, y = 2), genhz = "C")  # single-row group, must fall back internally
  )

  result <- estimate_correlation_matrix_robust(df, group_var = "genhz", min_group_n = 5)
  expect_equal(result$method, "empirical_grouped")
  expect_true(result$matrix["x", "y"] > 0.4)
  expect_true(all(eigen(result$matrix, only.values = TRUE)$values > 0))
})

test_that("estimate_correlation_matrix_robust edge case: n=1 total observation does not error", {
  df <- data.frame(x = 1, y = 2, z = 3)
  expect_no_error(result <- estimate_correlation_matrix_robust(df, min_group_n = 5))
  expect_equal(result$method, "global_fallback")
})

test_that("estimate_correlation_matrix_robust passing group_fallback_matrices=NULL/global_fallback_method='global_fallback' reproduces current behavior exactly (backward-compatibility guard)", {
  set.seed(19)
  n_per_group <- 30
  make_group <- function(rho) {
    x <- rnorm(n_per_group)
    y <- rho * x + rnorm(n_per_group, sd = sqrt(1 - rho^2))
    data.frame(x = x, y = y)
  }
  df <- rbind(
    cbind(make_group(0.7), genhz = "A"),
    cbind(make_group(0.7), genhz = "B"),
    cbind(data.frame(x = 1, y = 2), genhz = "C")
  )

  baseline <- estimate_correlation_matrix_robust(df, group_var = "genhz", min_group_n = 5)
  explicit_defaults <- estimate_correlation_matrix_robust(
    df, group_var = "genhz", min_group_n = 5,
    group_fallback_matrices = NULL, global_fallback_method = "global_fallback"
  )
  expect_identical(baseline, explicit_defaults)
})

test_that("estimate_correlation_matrix_robust substitutes a failing group's own fallback matrix instead of dropping it", {
  set.seed(23)
  n_per_group <- 30
  make_group <- function(rho) {
    x <- rnorm(n_per_group)
    y <- rho * x + rnorm(n_per_group, sd = sqrt(1 - rho^2))
    data.frame(x = x, y = y)
  }
  # Group "A" has enough data to estimate empirically (rho ~ 0.7); group "B"
  # has too little (n=2, below the default min_group_n=5) and has no
  # empirical result on its own.
  df <- rbind(
    cbind(make_group(0.7), genhz = "A"),
    cbind(data.frame(x = c(1, 2), y = c(1.1, 2.2)), genhz = "B")
  )
  fallback_b <- matrix(c(1, -0.5, -0.5, 1), nrow = 2, dimnames = list(c("x", "y"), c("x", "y")))

  result <- estimate_correlation_matrix_robust(
    df, group_var = "genhz", min_group_n = 5,
    group_fallback_matrices = list(B = fallback_b)
  )
  expect_equal(result$method, "empirical_grouped_kssl_blended")
  expect_true(all(eigen(result$matrix, only.values = TRUE)$values > 0))
  # Blended value should sit strictly between the two sources' correlations
  # (group A's empirical rho and group B's fallback -0.5), not equal either.
  expect_true(result$matrix["x", "y"] > -0.5 && result$matrix["x", "y"] < 0.9)
})

test_that("estimate_correlation_matrix_robust returns kssl_fallback_grouped when every group is fallback-sourced", {
  df <- rbind(
    cbind(data.frame(x = c(1, 2), y = c(1.1, 2.2)), genhz = "A"),
    cbind(data.frame(x = c(3, 4), y = c(2.9, 4.1)), genhz = "B")
  )
  fallback <- matrix(c(1, 0.4, 0.4, 1), nrow = 2, dimnames = list(c("x", "y"), c("x", "y")))

  result <- estimate_correlation_matrix_robust(
    df, group_var = "genhz", min_group_n = 5,
    group_fallback_matrices = list(A = fallback, B = fallback)
  )
  expect_equal(result$method, "kssl_fallback_grouped")
  expect_equal(result$matrix["x", "y"], 0.4)
})

test_that("estimate_correlation_matrix_robust still drops a group with no empirical data and no matching fallback entry", {
  set.seed(29)
  n_per_group <- 30
  make_group <- function(rho) {
    x <- rnorm(n_per_group)
    y <- rho * x + rnorm(n_per_group, sd = sqrt(1 - rho^2))
    data.frame(x = x, y = y)
  }
  df <- rbind(
    cbind(make_group(0.6), genhz = "A"),
    cbind(data.frame(x = c(1, 2), y = c(1.1, 2.2)), genhz = "B")  # too little data, no fallback entry for "B"
  )

  result <- estimate_correlation_matrix_robust(
    df, group_var = "genhz", min_group_n = 5,
    group_fallback_matrices = list(C = matrix(c(1, 0.9, 0.9, 1), nrow = 2))  # unrelated key, doesn't cover "B"
  )
  expect_equal(result$method, "empirical_grouped")
  expect_equal(unname(result$n_obs), n_per_group)  # only group A contributed
})

test_that("estimate_correlation_matrix_robust global_fallback_method overrides the final-tier method string", {
  df <- data.frame(x = 1, y = 2, z = 3)
  custom_fallback <- diag(3)
  rownames(custom_fallback) <- colnames(custom_fallback) <- c("x", "y", "z")

  result <- estimate_correlation_matrix_robust(
    df, min_group_n = 5, global_fallback = custom_fallback,
    global_fallback_method = "kssl_global_fallback"
  )
  expect_equal(result$method, "kssl_global_fallback")
  expect_equal(result$matrix, custom_fallback)
})

test_that("resolve_composition_groups activates the texture group only when all 3 members present", {
  config <- list(monte_carlo = list(composition_groups = list(
    texture = list(members = c("claytotal", "sandtotal", "silttotal"), pseudo = c("ilr1", "ilr2"))
  )))

  full <- resolve_composition_groups(c("dbovendry", "sandtotal", "claytotal", "silttotal"), config)
  expect_true(full$groups$texture$active)
  expect_equal(full$sim_properties, c("dbovendry", "ilr1", "ilr2"))

  partial <- resolve_composition_groups(c("dbovendry", "sandtotal", "claytotal"), config)
  expect_false(partial$groups$texture$active)
  expect_equal(partial$sim_properties, c("dbovendry", "sandtotal", "claytotal"))

  none <- resolve_composition_groups(c("dbovendry"), config)
  expect_false(none$groups$texture$active)
  expect_equal(none$sim_properties, c("dbovendry"))
})

test_that("restore_composition_properties inverse-transforms ilr1/ilr2 back to clay/sand/silt summing to 100", {
  config <- list(monte_carlo = list(composition_groups = list(
    texture = list(members = c("claytotal", "sandtotal", "silttotal"), pseudo = c("ilr1", "ilr2"))
  )))
  plan <- resolve_composition_groups(c("dbovendry", "sandtotal", "claytotal", "silttotal"), config)

  n_horizons <- 4; n_realizations <- 50
  set.seed(23)
  mu <- as.numeric(ilr_forward(20, 40, 40))
  sim <- array(
    NA_real_,
    dim = c(n_horizons, length(plan$sim_properties), n_realizations),
    dimnames = list(horizon = 1:n_horizons, property = plan$sim_properties, realization = 1:n_realizations)
  )
  sim[, "dbovendry", ] <- rnorm(n_horizons * n_realizations, 1.4, 0.05)
  for (h in seq_len(n_horizons)) {
    sim[h, "ilr1", ] <- rnorm(n_realizations, mu[1], 0.1)
    sim[h, "ilr2", ] <- rnorm(n_realizations, mu[2], 0.1)
  }

  restored <- restore_composition_properties(sim, plan$sim_properties, c("dbovendry", "sandtotal", "claytotal", "silttotal"), plan$groups)
  expect_equal(dimnames(restored)$property, c("dbovendry", "sandtotal", "claytotal", "silttotal"))
  expect_equal(restored[, "dbovendry", ], sim[, "dbovendry", ])

  totals <- restored[, "sandtotal", ] + restored[, "claytotal", ] + restored[, "silttotal", ]
  expect_equal(as.numeric(totals), rep(100, n_horizons * n_realizations), tolerance = 1e-6)
  expect_true(all(restored[, c("sandtotal", "claytotal", "silttotal"), ] >= 0))
  expect_true(all(restored[, c("sandtotal", "claytotal", "silttotal"), ] <= 100))
})
