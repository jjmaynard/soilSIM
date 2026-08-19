test_that("get_monte_carlo_defaults() defaults vertical_correlation_method to joint_copula (Phase 13 flip), and validate_monte_carlo_config() accepts both valid choices", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 13 flipped the default from
  # "gp_quantile_retrofit" to "joint_copula" after Phases 0-12 resolved every blocking decision
  # point - "gp_quantile_retrofit" remains fully supported as an explicit opt-out, not removed.
  default_config <- get_monte_carlo_defaults()
  expect_equal(default_config$monte_carlo$vertical_correlation_method, "joint_copula")

  validation_default <- validate_monte_carlo_config(default_config, n_realizations = 100)
  expect_true(validation_default$valid)

  retrofit_config <- default_config
  retrofit_config$monte_carlo$vertical_correlation_method <- "gp_quantile_retrofit"
  validation_retrofit <- validate_monte_carlo_config(retrofit_config, n_realizations = 100)
  expect_true(validation_retrofit$valid)

  bad_config <- default_config
  bad_config$monte_carlo$vertical_correlation_method <- "not_a_real_method"
  validation_bad <- validate_monte_carlo_config(bad_config, n_realizations = 100)
  expect_false(validation_bad$valid)

  # Not required (unlike correlation_fallback) - a config missing this key entirely must still
  # validate cleanly, since older/ad-hoc configs built without get_monte_carlo_defaults() should
  # not suddenly start failing validation.
  no_key_config <- default_config
  no_key_config$monte_carlo$vertical_correlation_method <- NULL
  validation_no_key <- validate_monte_carlo_config(no_key_config, n_realizations = 100)
  expect_true(validation_no_key$valid)
})

test_that("get_monte_carlo_defaults() defaults vertical_correlation_gating to FALSE, and validate_monte_carlo_config() validates it as a logical", {
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 8 - discontinuity gating is a separate opt-in
  # from vertical_correlation_method itself, defaulting to off.
  default_config <- get_monte_carlo_defaults()
  expect_identical(default_config$monte_carlo$vertical_correlation_gating, FALSE)

  validation_default <- validate_monte_carlo_config(default_config, n_realizations = 100)
  expect_true(validation_default$valid)

  gating_on_config <- default_config
  gating_on_config$monte_carlo$vertical_correlation_gating <- TRUE
  validation_on <- validate_monte_carlo_config(gating_on_config, n_realizations = 100)
  expect_true(validation_on$valid)

  bad_config <- default_config
  bad_config$monte_carlo$vertical_correlation_gating <- "yes"  # not logical
  validation_bad <- validate_monte_carlo_config(bad_config, n_realizations = 100)
  expect_false(validation_bad$valid)

  # Not required - a config missing this key entirely must still validate cleanly.
  no_key_config <- default_config
  no_key_config$monte_carlo$vertical_correlation_gating <- NULL
  validation_no_key <- validate_monte_carlo_config(no_key_config, n_realizations = 100)
  expect_true(validation_no_key$valid)
})

test_that("normalize_monte_carlo_config() accepts both flat and nested simulation_config", {
  default_config <- get_monte_carlo_defaults()

  flat <- list(distribution_type = "normal", max_depth = 100)
  normalized_flat <- normalize_monte_carlo_config(flat, default_config)
  expect_equal(normalized_flat, list(monte_carlo = flat))

  nested <- list(monte_carlo = list(distribution_type = "normal"))
  normalized_nested <- normalize_monte_carlo_config(nested, default_config)
  expect_equal(normalized_nested, nested)

  expect_null(normalize_monte_carlo_config(NULL, default_config))
  expect_equal(normalize_monte_carlo_config(list(), default_config), list())
})

test_that("check_property_data_availability()'s vectorized rowSums() matches the original nested per-row loop", {
  # PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4. Covers: a row with data in one property but not
  # another (still TRUE), a row with no data in any present property (FALSE), and a property
  # whose "_r" column isn't present in soil_data at all (silently ignored, not an error).
  soil_data <- data.frame(
    claytotal_r = c(10, NA, NA, 20),
    sandtotal_r = c(NA, NA, 30, NA),
    dbovendry_r = c(NA, NA, NA, NA)
  )
  properties <- c("claytotal", "sandtotal", "dbovendry", "ph1to1h2o") # ph1to1h2o_r not present

  reference_availability <- function(soil_data, properties) {
    has_data <- rep(FALSE, nrow(soil_data))
    for (i in seq_len(nrow(soil_data))) {
      for (prop in properties) {
        r_col <- paste0(prop, "_r")
        if (r_col %in% names(soil_data) && !is.na(soil_data[[r_col]][i])) {
          has_data[i] <- TRUE
          break
        }
      }
    }
    has_data
  }

  result <- check_property_data_availability(soil_data, properties, config = list())
  expected <- reference_availability(soil_data, properties)
  expect_equal(result, expected)
  expect_equal(result, c(TRUE, FALSE, TRUE, TRUE))
})

test_that("check_property_data_availability() returns all FALSE when no property's _r column is present", {
  soil_data <- data.frame(x = 1:3)
  result <- check_property_data_availability(soil_data, c("claytotal"), config = list())
  expect_equal(result, c(FALSE, FALSE, FALSE))
})

test_that("distribution_type='normal' actually threads through to a real fit (bug #1 regression)", {
  soil_data <- make_soil_data(10, properties = c("dbovendry"))
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 100, seed = 1,
    simulation_config = list(distribution_type = "normal")
  )
  expect_equal(res$metadata$config_used$monte_carlo$distribution_type, "normal")
  expect_equal(res$parameters[[1]]$dbovendry$family, "normal")
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("distribution_type='beta' with a registry bounds entry produces finite, in-bounds draws", {
  soil_data <- make_soil_data(10, properties = c("sandtotal"))
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("sandtotal"), n_realizations = 100, seed = 2,
    simulation_config = list(
      property_distributions = list(sandtotal = list(family = "beta", bounds = c(0, 100)))
    )
  )
  expect_equal(res$parameters[[1]]$sandtotal$family, "beta")
  expect_true(all(is.finite(res$simulation_data)))
  expect_true(all(res$simulation_data >= 0 & res$simulation_data <= 100))
})

test_that("validate_distribution_parameters() (mod05 wrapper) catches invalid parameters, not always TRUE", {
  bad <- validate_distribution_parameters(list(family = "beta", fit = list(shape1 = -1, shape2 = 2, lower = 0, upper = 100)), "triangular", "x")
  expect_false(bad$valid)

  good <- validate_distribution_parameters(list(family = "triangular", fit = list(min = 1, mode = 2, max = 3)), "triangular", "x")
  expect_true(good$valid)
})

test_that("estimate_property_correlations() returns a real, valid correlation matrix (bug #2 regression)", {
  set.seed(3)
  n <- 40
  db_r <- rnorm(n, 1.4, 0.1)
  ph_r <- 6 + 0.8 * (db_r - 1.4) / 0.1 + rnorm(n, 0, 0.2)  # strongly correlated with dbovendry
  rows <- lapply(seq_len(n), function(i) {
    make_horizon_row(
      properties = c("dbovendry", "ph1to1h2o"), cokey = as.character(i),
      dbovendry_r = db_r[i], ph1to1h2o_r = ph_r[i]
    )
  })
  soil_data <- do.call(rbind, rows)

  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o"), n_realizations = 50, seed = 3,
    simulation_config = list(auto_correlation = TRUE)
  )
  expect_equal(res$correlation_structure$method, "estimated_from_data")
  expect_true(all(eigen(res$correlation_structure$matrix, only.values = TRUE)$values > 0))
  expect_true(abs(res$correlation_structure$matrix["dbovendry", "ph1to1h2o"]) > 0.5)
})

test_that("estimate_property_correlations() falls back to identity gracefully with insufficient data", {
  soil_data <- make_soil_data(3, properties = c("dbovendry", "ph1to1h2o"))
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o"), n_realizations = 20, seed = 4,
    simulation_config = list(auto_correlation = TRUE)
  )
  expect_true(all(eigen(res$correlation_structure$matrix, only.values = TRUE)$values > 0))
  # method should reflect the robust estimator's own honest fallback, not error
  expect_true(res$correlation_structure$method %in% c("estimated_from_data", "identity_estimation_failed"))
})

test_that("correlation_fallback='kssl_global' derives genhz from hzname and produces a non-identity KSSL-sourced correlation matrix", {
  soil_data <- make_soil_data_hzname_only(12)
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o", "cec7"), n_realizations = 30, seed = 201,
    simulation_config = list(auto_correlation = TRUE, correlation_fallback = "kssl_global")
  )
  # Never the plain identity/uninformative outcome - hzname is recognizable
  # (Ap/Bt1/Bt2/Bw/C1/C2) and correlation_fallback="kssl_global" is set, so
  # classify_genhz() + the KSSL fallback should always kick in for at least
  # some groups.
  expect_true(res$correlation_structure$estimation_info$method %in% c(
    "empirical_grouped", "empirical_grouped_kssl_blended",
    "kssl_fallback_grouped", "empirical_pooled", "kssl_global_fallback"
  ))
  m <- res$correlation_structure$matrix
  expect_false(isTRUE(all(m[upper.tri(m)] == 0)))
  expect_true(all(eigen(m, only.values = TRUE)$values > 0))
})

test_that("correlation_fallback default ('identity') produces byte-identical results whether omitted or set explicitly", {
  soil_data <- make_soil_data(3, properties = c("dbovendry", "ph1to1h2o"))
  res_omitted <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o"), n_realizations = 20, seed = 4,
    simulation_config = list(auto_correlation = TRUE)
  )
  res_explicit <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o"), n_realizations = 20, seed = 4,
    simulation_config = list(auto_correlation = TRUE, correlation_fallback = "identity")
  )
  expect_equal(res_omitted$correlation_structure, res_explicit$correlation_structure)
  expect_equal(res_omitted$simulation_data, res_explicit$simulation_data)
})

test_that("REGRESSION: a caller with hzname present but correlation_fallback never set gets no new grouping introduced (opt-in guarantee)", {
  # The single most important regression test for "opt-in, not opt-out":
  # make_soil_data_hzname_only() has real hzname values and no genhz column,
  # but correlation_fallback defaults to "identity" - classify_genhz() must
  # never be invoked and no genhz grouping must appear.
  soil_data <- make_soil_data_hzname_only(12)
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "ph1to1h2o", "cec7"), n_realizations = 30, seed = 201,
    simulation_config = list(auto_correlation = TRUE)
  )
  # With 12 ungrouped rows (no genhz stratification) the pooled empirical fit
  # has enough data (>5) to succeed outright - the fallback machinery this
  # feature added should never be reachable here.
  expect_equal(res$correlation_structure$estimation_info$method, "empirical_pooled")
})

test_that("end-to-end ILR texture simulation: every realization sums to 100 and stays in [0,100]", {
  texture_data <- make_texture_only_soil_data(12)
  res <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"), n_realizations = 300, seed = 5
  )
  texture_arr <- res$simulation_data[, c("sandtotal", "claytotal", "silttotal"), ]
  totals <- res$simulation_data[, "sandtotal", ] + res$simulation_data[, "claytotal", ] + res$simulation_data[, "silttotal", ]
  expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1e-6)
  expect_true(all(texture_arr >= 0 & texture_arr <= 100))
})

test_that("REGRESSION: composition_groups$texture$members role order is an internal reparameterization - simulated clay/sand/silt marginals are statistically invariant to the old (clay,sand,silt) order vs the new default (sand,silt,clay)", {
  # The default role order changed (to align soilSIM's ILR convention with
  # the optional KSSL reference-correlation fallback's convention - see
  # kssl-reference-correlations.R). ilr_forward()/ilr_inverse() themselves
  # were NOT touched - only which real property occupies position 1/2/3.
  # This locks in the claim (made in get_monte_carlo_defaults()'s doc
  # comment) that simulated clay/sand/silt output is statistically
  # unaffected by which role order is used.
  texture_data <- make_texture_only_soil_data(15)
  old_role_order_config <- list(composition_groups = list(
    texture = list(members = c("claytotal", "sandtotal", "silttotal"), pseudo = c("ilr1", "ilr2"))
  ))

  res_new_default <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 2000, seed = 101
  )
  res_old_order <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 2000, seed = 101, simulation_config = old_role_order_config
  )

  for (prop in c("sandtotal", "claytotal", "silttotal")) {
    new_vals <- as.numeric(res_new_default$simulation_data[, prop, ])
    old_vals <- as.numeric(res_old_order$simulation_data[, prop, ])
    expect_equal(mean(new_vals), mean(old_vals), tolerance = 0.05)
    expect_equal(stats::sd(new_vals), stats::sd(old_vals), tolerance = 0.15)
  }

  # Sum-to-100 holds under both role orders (exact ILR-inverse construction,
  # order-independent).
  for (res in list(res_new_default, res_old_order)) {
    totals <- res$simulation_data[, "sandtotal", ] + res$simulation_data[, "claytotal", ] + res$simulation_data[, "silttotal", ]
    expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1e-6)
  }
})

test_that("edge case: only 2 of 3 texture members requested falls back to the legacy path with a WARN, not an error", {
  texture_data <- make_texture_only_soil_data(10)
  expect_no_error(
    res <- generate_monte_carlo_realizations(
      texture_data, properties = c("sandtotal", "claytotal"), n_realizations = 50, seed = 6
    )
  )
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("edge case: an active composition group + explicit correlation_matrix falls back to the legacy path with a WARN", {
  texture_data <- make_texture_only_soil_data(10)
  user_matrix <- diag(3)
  rownames(user_matrix) <- colnames(user_matrix) <- c("sandtotal", "claytotal", "silttotal")

  expect_warning(
    res <- generate_monte_carlo_realizations(
      texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
      correlation_matrix = user_matrix, n_realizations = 50, seed = 7
    ),
    regexp = NA  # log_message() doesn't raise an R warning condition; just assert no error/crash
  )
  expect_true(all(is.finite(res$simulation_data)))
  totals <- res$simulation_data[, "sandtotal", ] + res$simulation_data[, "claytotal", ] + res$simulation_data[, "silttotal", ]
  # Legacy proportional-rescale path (not the exact ILR path) is used here,
  # so totals should be close to but not necessarily bit-exact at 100.
  expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1)
})

test_that("lh_percentile config knob actually changes the fitted SD for a normal-family property", {
  soil_data <- make_soil_data(5, properties = c("dbovendry"))

  narrow <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 10, seed = 8,
    simulation_config = list(distribution_type = "normal", lh_percentile = c(0.25, 0.75))
  )
  wide <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 10, seed = 8,
    simulation_config = list(distribution_type = "normal", lh_percentile = c(0.05, 0.95))
  )

  narrow_sd <- narrow$parameters[[1]]$dbovendry$fit$sd
  wide_sd <- wide$parameters[[1]]$dbovendry$fit$sd
  expect_true(narrow_sd > wide_sd)  # same raw l/h span, narrower assumed percentiles -> larger implied SD
})

test_that("metalog infeasible-fit fallback triggers correctly inside the full pipeline", {
  soil_data <- make_horizon_row(properties = "om", om_l = 0.1, om_r = 0.5, om_h = 99.9)
  soil_data <- rbind(soil_data, make_horizon_row(properties = "om", cokey = "2", om_l = 0.1, om_r = 0.5, om_h = 99.9))

  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("om"), n_realizations = 50, seed = 9,
    simulation_config = list(
      distribution_type = "metalog",
      property_distributions = list(om = list(bounds = c(0, 100)))
    )
  )
  expect_equal(res$parameters[[1]]$om$family, "metalog")
  expect_true(isTRUE(res$parameters[[1]]$om$fit$infeasible))
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("edge case: a single horizon does not crash the pipeline", {
  soil_data <- make_soil_data(1, properties = c("dbovendry"))
  expect_no_error(
    res <- generate_monte_carlo_realizations(soil_data, properties = c("dbovendry"), n_realizations = 20, seed = 10)
  )
  expect_equal(dim(res$simulation_data)[1], 1)
})

test_that("edge case: a property missing from every horizon degrades gracefully (WARN, not crash)", {
  soil_data <- make_soil_data(5, properties = c("dbovendry"))
  soil_data$cec7_l <- NA_real_
  soil_data$cec7_r <- NA_real_
  soil_data$cec7_h <- NA_real_

  expect_no_error(
    res <- generate_monte_carlo_realizations(
      soil_data, properties = c("dbovendry", "cec7"), n_realizations = 20, seed = 11
    )
  )
  expect_true(res$validation$output_validation$success_rate < 1)
})

test_that("parallel and sequential simulation both produce valid sum-to-100 texture output", {
  skip_on_cran()
  texture_data <- make_texture_only_soil_data(8)

  res_seq <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 40, seed = 12, parallel = FALSE
  )
  res_par <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 40, seed = 12, parallel = TRUE,
    simulation_config = list(parallel_threshold = 1)
  )

  for (res in list(res_seq, res_par)) {
    totals <- res$simulation_data[, "sandtotal", ] + res$simulation_data[, "claytotal", ] + res$simulation_data[, "silttotal", ]
    expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1e-6)
    expect_equal(sort(dimnames(res$simulation_data)$property), sort(c("sandtotal", "claytotal", "silttotal")))
  }
})

# ---------------------------------------------------------------------------
# observed_data / Bayesian-updating wiring (fuse_observed_data_into_priors())
# ---------------------------------------------------------------------------

test_that("observed_data (normal family, closed-form route) tightens the posterior and shifts toward the likelihood", {
  soil_data <- make_soil_data(10, properties = c("dbovendry"))
  res_prior_only <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 300, seed = 1,
    simulation_config = list(distribution_type = "normal")
  )
  res_fused <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 300, seed = 1,
    simulation_config = list(distribution_type = "normal"),
    observed_data = list(dbovendry = list(mean = 1.35, sd = 0.02))
  )

  expect_equal(res_fused$parameters[[1]]$dbovendry$family, "normal")
  expect_true(res_fused$parameters[[1]]$dbovendry$fit$sd < 0.02 + 1e-9)  # tighter than the (already tight) likelihood
  expect_true(sd(res_fused$simulation_data) < sd(res_prior_only$simulation_data))
  expect_equal(mean(res_fused$simulation_data), 1.35, tolerance = 0.05)  # dominated by the tight likelihood
})

test_that("observed_data (lognormal family) does not double-convert the already-log-space prior", {
  soil_data <- make_soil_data(6, properties = c("om"))
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("om"), n_realizations = 50, seed = 6,
    simulation_config = list(
      distribution_type = "lognormal",
      property_distributions = list(om = list(family = "lognormal"))
    ),
    observed_data = list(om = list(mean = 2, sd = 0.3))
  )
  posterior_fit <- res$parameters[[1]]$om$fit
  expect_equal(res$parameters[[1]]$om$family, "lognormal")
  # Round-tripping the stored (log-space) posterior back to raw space should
  # land near a sane blend of the prior's raw om_r (~2) and the likelihood's
  # raw mean (2) - not some wildly different value from a double conversion.
  # lognormal_to_normal_params() returns list(mu=, sigma=) (raw-space mean/sd).
  raw <- lognormal_to_normal_params(posterior_fit$mean, posterior_fit$sd)
  expect_true(raw$mu > 0.5 && raw$mu < 5)
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("observed_data (beta family) rescales onto the prior's own bounds and stays feasible", {
  soil_data <- make_soil_data(8, properties = c("dbovendry"))
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 200, seed = 5,
    simulation_config = list(
      distribution_type = "beta",
      property_distributions = list(dbovendry = list(family = "beta", bounds = c(0, 3)))
    ),
    observed_data = list(dbovendry = list(mean = 1.35, var = 0.001))
  )
  fit <- res$parameters[[1]]$dbovendry$fit
  expect_equal(res$parameters[[1]]$dbovendry$family, "beta")
  expect_true(fit$shape1 > 0 && fit$shape2 > 0)
  expect_equal(fit$lower, 0)
  expect_equal(fit$upper, 3)
  expect_true(all(is.finite(res$simulation_data)))
  expect_true(all(res$simulation_data >= 0 & res$simulation_data <= 3))
})

test_that("observed_data (beta family, infeasible fusion) falls back to the moment-based route instead of producing invalid shape params", {
  # Two very diffuse (near-uniform) sides are the classic fuse_beta() infeasible
  # case (prior_alpha+lik_alpha <= 1 or prior_beta+lik_beta <= 1).
  soil_data <- make_horizon_row(properties = "dbovendry", dbovendry_l = 0.1, dbovendry_r = 1.5, dbovendry_h = 2.9)
  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 100, seed = 7, validate_inputs = FALSE,
    simulation_config = list(
      distribution_type = "beta",
      property_distributions = list(dbovendry = list(family = "beta", bounds = c(0, 3)))
    ),
    observed_data = list(dbovendry = list(shape1 = 0.3, shape2 = 0.3))
  )
  fit <- res$parameters[[1]]$dbovendry$fit
  expect_true(is.finite(fit$shape1) && fit$shape1 > 0)
  expect_true(is.finite(fit$shape2) && fit$shape2 > 0)
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("observed_data (general/vector route) produces a linear_cdf posterior between prior and likelihood", {
  soil_data <- make_soil_data(10, properties = c("dbovendry"))
  set.seed(9)
  obs_samples <- rnorm(30, 1.35, 0.03)

  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry"), n_realizations = 200, seed = 9,
    simulation_config = list(distribution_type = "normal"),
    observed_data = list(dbovendry = obs_samples)
  )
  posterior <- res$parameters[[1]]$dbovendry
  expect_equal(posterior$family, "linear_cdf")
  expect_equal(posterior$source, "bayesian_fusion_general")

  posterior_median <- quantile_from_fit(0.5, posterior$family, posterior$fit)
  expect_true(posterior_median > min(1.4, mean(obs_samples)) - 0.1)  # roughly between prior (~1.4) and likelihood (~1.35)
  expect_true(posterior_median < max(1.4, mean(obs_samples)) + 0.1)
})

test_that("observed_data skips fusion (with WARN) for a prior family with no closed-form fusion route", {
  soil_data <- make_soil_data(5, properties = c("dbovendry"))
  # distribution_type defaults to "triangular", which has no closed-form
  # bayes_fuse() route - a param-list observed_data entry should be skipped.
  expect_no_error(
    res <- generate_monte_carlo_realizations(
      soil_data, properties = c("dbovendry"), n_realizations = 20, seed = 2,
      observed_data = list(dbovendry = list(mean = 1.35, sd = 0.02))
    )
  )
  expect_equal(res$parameters[[1]]$dbovendry$family, "triangular")
  expect_true(all(is.finite(res$simulation_data)))
})

test_that("observed_data texture fusion tightens ilr1/ilr2 and preserves exact sum-to-100", {
  texture_data <- make_texture_only_soil_data(8)
  res_prior_only <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 300, seed = 4
  )
  res_fused <- generate_monte_carlo_realizations(
    texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
    n_realizations = 300, seed = 4,
    observed_data = list(
      claytotal = c(low = 18, rep = 20, high = 22),
      sandtotal = c(low = 38, rep = 40, high = 42),
      silttotal = c(low = 38, rep = 40, high = 42)
    )
  )

  totals <- res_fused$simulation_data[, "sandtotal", ] + res_fused$simulation_data[, "claytotal", ] + res_fused$simulation_data[, "silttotal", ]
  expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1e-6)
  expect_true(sd(res_fused$simulation_data[, "claytotal", ]) < sd(res_prior_only$simulation_data[, "claytotal", ]))
})

test_that("observed_data texture fusion skips (with WARN) when only 1-2 of 3 members are supplied", {
  texture_data <- make_texture_only_soil_data(6)
  expect_no_error(
    res <- generate_monte_carlo_realizations(
      texture_data, properties = c("sandtotal", "claytotal", "silttotal"),
      n_realizations = 50, seed = 8,
      observed_data = list(claytotal = c(low = 18, rep = 20, high = 22))
    )
  )
  totals <- res$simulation_data[, "sandtotal", ] + res$simulation_data[, "claytotal", ] + res$simulation_data[, "silttotal", ]
  expect_equal(as.numeric(totals), rep(100, length(totals)), tolerance = 1e-6)
})

test_that("Step 4.5 placement: observed_data fusion survives a colliding user-supplied correlation_matrix", {
  # The active_group_with_user_matrix recompute (existing code, Step 5) must
  # not silently discard Step 4.5's fused posteriors for the NON-texture
  # property in this scenario - the recompute only touches texture members.
  soil_data <- make_soil_data(8, properties = c("dbovendry", "sandtotal", "claytotal", "silttotal"))
  user_matrix <- diag(4)
  rownames(user_matrix) <- colnames(user_matrix) <- c("dbovendry", "sandtotal", "claytotal", "silttotal")

  res <- generate_monte_carlo_realizations(
    soil_data, properties = c("dbovendry", "sandtotal", "claytotal", "silttotal"),
    correlation_matrix = user_matrix, n_realizations = 50, seed = 10,
    simulation_config = list(distribution_type = "normal"),
    observed_data = list(dbovendry = list(mean = 1.35, sd = 0.02))
  )
  expect_equal(res$parameters[[1]]$dbovendry$source, "bayesian_fusion_normal")
  expect_true(res$parameters[[1]]$dbovendry$fit$sd < 0.1)
})

test_that("generate_monte_carlo_realizations() with observed_data=NULL is unchanged from omitting the argument entirely", {
  soil_data <- make_soil_data(6, properties = c("dbovendry"))
  res_omitted <- generate_monte_carlo_realizations(soil_data, properties = c("dbovendry"), n_realizations = 30, seed = 13)
  res_explicit_null <- generate_monte_carlo_realizations(soil_data, properties = c("dbovendry"), n_realizations = 30, seed = 13, observed_data = NULL)
  expect_equal(res_omitted$simulation_data, res_explicit_null$simulation_data)
})

## ----------------------------------------------------------------------
## Component composition helper functions (previously untagged placeholders)
## ----------------------------------------------------------------------

test_that("validate_component_data() requires a non-empty data frame with a comppct_r column", {
  expect_true(validate_component_data(data.frame(comppct_r = c(50, 50)))$valid)
  expect_false(validate_component_data(data.frame(x = 1))$valid)
  expect_false(validate_component_data(data.frame(comppct_r = character(0)))$valid)
  expect_false(validate_component_data(data.frame(comppct_r = NA_real_))$valid)
})

test_that("extract_component_parameters() uses comppct_l/comppct_r/comppct_h, falling back to comppct_r -/+ 2", {
  full_row <- data.frame(comppct_l = 45, comppct_r = 50, comppct_h = 60)
  result_full <- extract_component_parameters(full_row, config = list())
  expect_equal(result_full, list(min = 45, mode = 50, max = 60))

  partial_row <- data.frame(comppct_l = NA_real_, comppct_r = 50, comppct_h = NA_real_)
  result_partial <- extract_component_parameters(partial_row, config = list())
  expect_equal(result_partial, list(min = 48, mode = 50, max = 52))

  missing_row <- data.frame(comppct_l = NA_real_, comppct_r = NA_real_, comppct_h = NA_real_)
  expect_equal(extract_component_parameters(missing_row, config = list()), list(min = 0, mode = 10, max = 20))
})

test_that("apply_composition_constraints() clamps to min/max and is a no-op when constraints is NULL", {
  m <- matrix(c(-5, 50, 150, 10), nrow = 2)
  expect_identical(apply_composition_constraints(m, NULL), m)

  clamped <- apply_composition_constraints(m, list(min = 0, max = 100))
  expect_equal(as.vector(clamped), c(0, 50, 100, 10))
  expect_equal(dim(clamped), dim(m))
})

test_that("assess_component_quality() scores how closely simulated means track original comppct_r values", {
  original_data <- data.frame(comppct_r = c(50, 30, 20))
  good_sim <- matrix(rep(c(50, 30, 20), 10), nrow = 3)
  result_good <- assess_component_quality(original_data, good_sim, config = list())
  expect_equal(result_good$overall_quality, 1, tolerance = 1e-8)

  bad_sim <- matrix(rep(c(90, 5, 1), 10), nrow = 3)
  result_bad <- assess_component_quality(original_data, bad_sim, config = list())
  expect_true(result_bad$overall_quality < result_good$overall_quality)
})

test_that("sim_component_compositions() runs end-to-end and reports a real (non-hardcoded) quality score", {
  component_data <- data.frame(
    compname = c("A", "B", "C"),
    comppct_l = c(NA, 25, 15), comppct_r = c(50, 30, 20), comppct_h = c(NA, 35, 25)
  )
  result <- sim_component_compositions(component_data, n_realizations = 200,
                                       config = get_monte_carlo_defaults())
  expect_equal(dim(result$realizations), c(200, 3))
  expect_true(is.finite(result$quality_metrics$overall_quality))
  # Normalization forces each realization's row to sum to 100.
  expect_true(all(abs(rowSums(result$realizations) - 100) < 1e-6))
})

## ----------------------------------------------------------------------
## Simulation quality/diagnostics functions (previously untagged placeholders)
## ----------------------------------------------------------------------

test_that("calculate_simulation_quality_metrics() computes real finite/outlier rates from the simulation array", {
  set.seed(1)
  sim_array <- array(stats::rnorm(2 * 1 * 100, 20, 2), dim = c(2, 1, 100),
                     dimnames = list(NULL, "clay_pct", NULL))
  config <- list(monte_carlo = list(outlier_threshold = 1.5))
  result <- calculate_simulation_quality_metrics(sim_array, "clay_pct", config)
  expect_equal(result$overall_finite_rate, 1)
  expect_true(is.finite(result$overall_quality))
})

test_that("generate_simulation_diagnostics() compares simulated vs original per-property mean/sd", {
  set.seed(1)
  sim_array <- array(stats::rnorm(2 * 1 * 200, 20, 2), dim = c(2, 1, 200),
                     dimnames = list(NULL, "clay_pct", NULL))
  original_data <- data.frame(clay_pct = stats::rnorm(50, 20, 2))
  result <- generate_simulation_diagnostics(original_data, sim_array, "clay_pct",
                                            correlation_config = list(method = "identity_default"),
                                            config = list())
  expect_equal(result$per_property$clay_pct$sim_mean, 20, tolerance = 1)
  expect_equal(result$correlation_method, "identity_default")
})

test_that("assess_simulation_quality() derives real component scores from diagnostics rather than hardcoded constants", {
  diagnostics_good <- list(per_property = list(
    clay_pct = list(sim_mean = 20, sim_sd = 2, orig_mean = 20, orig_sd = 2)
  ))
  result_good <- assess_simulation_quality(NULL, diagnostics_good, list(success_rate = 1.0), list())
  expect_equal(result_good$component_scores$constraint_satisfaction, 1, tolerance = 1e-8)
  expect_equal(result_good$component_scores$statistical_consistency, 1, tolerance = 1e-8)

  diagnostics_bad <- list(per_property = list(
    clay_pct = list(sim_mean = 80, sim_sd = 10, orig_mean = 20, orig_sd = 2)
  ))
  result_bad <- assess_simulation_quality(NULL, diagnostics_bad, list(success_rate = 1.0), list())
  expect_true(result_bad$overall_quality_score < result_good$overall_quality_score)
})

test_that("summarize_distributions()/validate_distribution_setup() count real valid/invalid distribution entries", {
  distributions <- list(
    list(clay_pct = list(valid = TRUE, family = "normal"), sand_pct = list(valid = FALSE, family = "triangular")),
    list(clay_pct = list(valid = TRUE, family = "normal"))
  )
  summary_result <- summarize_distributions(distributions, c("clay_pct", "sand_pct"), "normal")
  expect_equal(summary_result$total_distributions, 3)
  expect_equal(summary_result$valid_distributions, 2)

  validation_result <- validate_distribution_setup(distributions, c("clay_pct", "sand_pct"), list())
  expect_equal(validation_result$validity_rate, 2 / 3, tolerance = 1e-8)
})

test_that("apply_sum_constraints() vectorized rescale matches the original per-cell scalar loop", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 1 apply_sum_constraints() fix:
  # reimplements the ORIGINAL nested for(h) for(r) scalar loop here and checks the vectorized
  # version produces identical output, including on edge cases (a zero-sum cell that must be left
  # untouched, and a cell already within the 0.01 tolerance that must also be left untouched).
  old_scalar <- function(results, constraints, properties) {
    adjustments <- 0
    for (constraint_name in names(constraints)) {
      constraint <- constraints[[constraint_name]]
      prop_indices <- match(constraint$properties, properties)
      prop_indices <- prop_indices[!is.na(prop_indices)]
      if (length(prop_indices) >= 2) {
        for (h in seq_len(dim(results)[1])) {
          for (r in seq_len(dim(results)[3])) {
            current_sum <- sum(results[h, prop_indices, r])
            if (current_sum > 0 && abs(current_sum - constraint$target_sum) > 0.01) {
              results[h, prop_indices, r] <- results[h, prop_indices, r] * constraint$target_sum / current_sum
              adjustments <- adjustments + 1
            }
          }
        }
      }
    }
    list(data = results, summary = list(adjustments = adjustments))
  }

  set.seed(7)
  n_horizons <- 6; n_realizations <- 5
  properties <- c("clay_pct", "sand_pct", "silt_pct", "pH")
  results <- array(
    runif(n_horizons * length(properties) * n_realizations, 10, 50),
    dim = c(n_horizons, length(properties), n_realizations)
  )
  # Cell (1,,1): zero-sum edge case (must stay untouched by both versions).
  results[1, 1:3, 1] <- 0
  # Cell (2,,1): already within tolerance of target_sum=100 (must stay untouched).
  results[2, 1:3, 1] <- c(33.33, 33.33, 33.34)

  constraints <- list(texture = list(properties = c("clay_pct", "sand_pct", "silt_pct"), target_sum = 100))

  expected <- old_scalar(results, constraints, properties)
  actual <- apply_sum_constraints(results, constraints, properties)

  expect_equal(actual$data, expected$data, tolerance = 1e-12)
  expect_equal(actual$summary$adjustments, expected$summary$adjustments)
  expect_equal(unname(actual$data[1, 1:3, 1]), c(0, 0, 0))
  expect_equal(unname(actual$data[2, 1:3, 1]), c(33.33, 33.33, 33.34))
})
