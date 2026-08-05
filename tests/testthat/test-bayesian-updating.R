test_that("bayes_update_normal_normal() satisfies its documented invariants", {
  post <- bayes_update_normal_normal(prior_mu = 10, prior_sigma = 2, lik_mu = 12, lik_sigma = 1)
  expect_true(post$sigma <= min(2, 1))
  expect_true(post$mu >= min(10, 12) && post$mu <= max(10, 12))

  swapped <- bayes_update_normal_normal(prior_mu = 12, prior_sigma = 1, lik_mu = 10, lik_sigma = 2)
  expect_equal(post$mu, swapped$mu, tolerance = 1e-9)
  expect_equal(post$sigma, swapped$sigma, tolerance = 1e-9)
})

test_that("bayes_update_normal_normal() posterior sigma is strictly tighter when informative on both sides", {
  post <- bayes_update_normal_normal(10, 2, 12, 1)
  expect_true(post$sigma < 1)
})

test_that("normal_to_lognormal_params()/lognormal_to_normal_params() round trip", {
  mu <- 5; sigma <- 1.2
  log_params <- normal_to_lognormal_params(mu, sigma)
  back <- lognormal_to_normal_params(log_params$mu, log_params$sigma)
  expect_equal(back$mu, mu, tolerance = 1e-9)
  expect_equal(back$sigma, sigma, tolerance = 1e-9)
})

test_that("fuse_beta()/fuse_gamma() feasibility flags correctly identify infeasible combinations", {
  feasible <- fuse_beta(prior_alpha = 3, prior_beta = 5, lik_alpha = 2, lik_beta = 4)
  expect_true(feasible$feasible)
  infeasible <- fuse_beta(prior_alpha = 0.3, prior_beta = 0.3, lik_alpha = 0.3, lik_beta = 0.3)
  expect_false(infeasible$feasible)

  feasible_g <- fuse_gamma(prior_shape = 3, prior_rate = 1, lik_shape = 2, lik_rate = 1)
  expect_true(feasible_g$feasible)
  infeasible_g <- fuse_gamma(prior_shape = 0.2, prior_rate = 1, lik_shape = 0.2, lik_rate = 1)
  expect_false(infeasible_g$feasible)
})

test_that("fuse_beta() conjugate result matches grid-based numerical Bayesian updating", {
  set.seed(21)
  prior_samples <- rbeta(20000, 5, 3)
  lik_samples <- rbeta(20000, 4, 6)
  fused <- fuse_beta(5, 3, 4, 6)
  expect_true(fused$feasible)

  numeric_post <- bayesian_update(prior_samples, lik_samples, grid_range = c(0, 1), grid_resolution = 0.001)
  expect_equal(mean(numeric_post), fused$alpha / (fused$alpha + fused$beta), tolerance = 0.02)
})

test_that("moments_to_beta()/beta_to_moments() and moments_to_gamma()/gamma_to_moments() round trip", {
  beta_params <- moments_to_beta(mean = 0.3, var = 0.02)
  back_moments <- beta_to_moments(beta_params$alpha, beta_params$beta)
  expect_equal(back_moments$mean, 0.3, tolerance = 1e-6)
  expect_equal(back_moments$var, 0.02, tolerance = 1e-6)

  gamma_params <- moments_to_gamma(mean = 5, var = 2)
  back_gamma <- gamma_to_moments(gamma_params$shape, gamma_params$rate)
  expect_equal(back_gamma$mean, 5, tolerance = 1e-6)
  expect_equal(back_gamma$var, 2, tolerance = 1e-6)
})

test_that("bayes_fuse() dispatches correctly by family", {
  n <- bayes_fuse(list(mu = 10, sigma = 2), list(mu = 12, sigma = 1), family = "normal")
  expect_equal(n, bayes_update_normal_normal(10, 2, 12, 1))

  b <- bayes_fuse(list(alpha = 3, beta = 5), list(alpha = 2, beta = 4), family = "beta")
  expect_equal(b, fuse_beta(3, 5, 2, 4))

  g <- bayes_fuse(list(shape = 3, rate = 1), list(shape = 2, rate = 1), family = "gamma")
  expect_equal(g, fuse_gamma(3, 1, 2, 1))
})

test_that("bayesian_update() general route posterior mean lands between prior and likelihood means", {
  set.seed(23)
  prior_samples <- rnorm(2000, 5, 1)
  lik_samples <- rnorm(2000, 8, 1)
  post <- bayesian_update(prior_samples, lik_samples, n = 2000)
  expect_true(mean(post) > 5 && mean(post) < 8)
})

test_that("bayesian_update()'s general route roughly reproduces bayes_update_normal_normal() for two actual Normals", {
  set.seed(29)
  prior_samples <- rnorm(20000, 10, 2)
  lik_samples <- rnorm(20000, 12, 1)
  closed_form <- bayes_update_normal_normal(10, 2, 12, 1)
  general <- bayesian_update(prior_samples, lik_samples, grid_resolution = 0.02, n = 5000)
  expect_equal(mean(general), closed_form$mu, tolerance = 0.15)
  expect_equal(sd(general), closed_form$sigma, tolerance = 0.15)
})

test_that("bayesian_update() respects its n argument", {
  set.seed(31)
  post <- bayesian_update(rnorm(500), rnorm(500), n = 250)
  expect_length(post, 250)
})

test_that("fuse_bivariate_normal() matches a Monte Carlo density-multiplication cross-check", {
  mu1 <- c(0, 0); Sigma1 <- diag(c(1, 1))
  mu2 <- c(1, 1); Sigma2 <- diag(c(0.5, 0.5))
  fused <- fuse_bivariate_normal(mu1, Sigma1, mu2, Sigma2)

  # Cross-check: precision-weighted mean should sit closer to the more
  # confident (smaller-variance) side (mu2), and posterior variance should
  # be smaller than both inputs'.
  expect_true(fused$mu[1] > mu1[1] && fused$mu[1] < mu2[1] * 1.01)
  expect_true(all(diag(fused$Sigma) < diag(Sigma1)))
  expect_true(all(diag(fused$Sigma) < diag(Sigma2)))
})

test_that("fuse_texture_group_from_triplets() produces posterior samples summing to 100 and staying in bounds", {
  set.seed(37)
  prior_triplets <- list(clay = c(15, 20, 25), sand = c(35, 40, 45), silt = c(35, 40, 45))
  lik_triplets <- list(clay = c(20, 25, 30), sand = c(30, 35, 40), silt = c(35, 40, 45))

  result <- fuse_texture_group_from_triplets(
    prior_triplets, lik_triplets, z_prior = qnorm(0.95), z_lik = qnorm(0.95), n_samples = 500
  )
  expect_equal(rowSums(result$posterior_samples), rep(100, 500), tolerance = 1e-6)
  expect_true(all(result$posterior_samples >= 0 & result$posterior_samples <= 100))
})

test_that("fuse_texture_group_from_triplets() avoids the sum-to-100 violation independent per-fraction fusion has", {
  set.seed(41)
  prior_triplets <- list(clay = c(10, 15, 20), sand = c(50, 60, 70), silt = c(20, 25, 30))
  lik_triplets <- list(clay = c(15, 20, 25), sand = c(40, 50, 60), silt = c(25, 30, 35))

  joint <- fuse_texture_group_from_triplets(prior_triplets, lik_triplets, qnorm(0.95), qnorm(0.95), n_samples = 200)
  expect_equal(rowSums(joint$posterior_samples), rep(100, 200), tolerance = 1e-6)

  # Independent per-fraction fusion (naive baseline): fit a Beta to each side/
  # fraction, fuse independently, and show the resulting point estimates do
  # NOT sum to 100 - the exact defect the joint ILR path exists to avoid.
  fit_beta_from_triplet <- function(triplet) {
    fit_beta_mle_newton(triplet, bounds = c(0, 100))
  }
  independent_alpha_beta <- lapply(list(clay = list(prior_triplets$clay, lik_triplets$clay),
                                        sand = list(prior_triplets$sand, lik_triplets$sand),
                                        silt = list(prior_triplets$silt, lik_triplets$silt)),
                                    function(pair) {
                                      prior_fit <- fit_beta_from_triplet(pair[[1]])
                                      lik_fit <- fit_beta_from_triplet(pair[[2]])
                                      fuse_beta(prior_fit$shape1, prior_fit$shape2, lik_fit$shape1, lik_fit$shape2)
                                    })
  independent_means <- vapply(independent_alpha_beta, function(f) f$alpha / (f$alpha + f$beta), numeric(1)) * 100
  expect_false(isTRUE(all.equal(sum(independent_means), 100, tolerance = 1e-6)))
})

test_that("fuse_property() dispatches on input shape and errors on mismatched method override", {
  set.seed(43)
  general_result <- fuse_property(rnorm(200, 5, 1), rnorm(200, 8, 1))
  expect_true(is.numeric(general_result))

  closed_form_result <- fuse_property(list(mu = 10, sigma = 2), list(mu = 12, sigma = 1), family = "normal")
  expect_equal(closed_form_result, bayes_update_normal_normal(10, 2, 12, 1))

  expect_error(fuse_property(rnorm(10), list(mu = 1, sigma = 1)), "SAME shape")
  expect_error(fuse_property(list(mu = 1, sigma = 1), list(mu = 2, sigma = 1), family = "normal", method = "general"), "implies")
  expect_error(fuse_property(list(alpha = 1), list(alpha = 2)), "family is required")
})

test_that("structural boundary: bayesian-updating.R never calls back into monte-carlo.R", {
  # bayesian-updating.R must stay genuinely standalone/independently testable
  # - its own primitives never change to accommodate the pipeline wiring, the
  # pipeline (monte-carlo.R) is the only side that knows about the bridge.
  # Reading source .R files by relative path isn't reliable once the package
  # is installed (R CMD check runs tests against the installed/lazy-loaded
  # package, not the source tree), so this inspects the loaded function
  # BODIES via deparse(body()), which works identically under
  # devtools::load_all() or a real install.
  bayesian_updating_exports <- c(
    "bayes_update_normal_normal", "normal_to_lognormal_params", "lognormal_to_normal_params",
    "fuse_beta", "fuse_gamma", "moments_to_beta", "moments_to_gamma", "beta_to_moments",
    "gamma_to_moments", "bayes_fuse", "bayesian_update", "fuse_bivariate_normal",
    "fuse_texture_group_from_triplets", "fuse_property"
  )
  monte_carlo_only_functions <- c(
    "generate_monte_carlo_realizations", "simulate_correlated_properties", "sim_component_compositions",
    "setup_distributions", "prepare_simulation_parameters", "configure_correlation_structure",
    "apply_simulation_constraints", "extract_property_parameters", "estimate_property_correlations",
    "prepare_simulation_data", "run_sequential_simulation", "run_parallel_simulation",
    "validate_monte_carlo_inputs", "validate_simulation_output", "get_monte_carlo_defaults",
    "normalize_monte_carlo_config", "get_constraint_rules", "get_sum_constraints",
    "as_lrh_triplet"
  )

  for (bu_fn in bayesian_updating_exports) {
    fn_obj <- get(bu_fn, envir = asNamespace("soilSIM"))
    fn_src <- paste(deparse(body(fn_obj)), collapse = "\n")
    for (mc_fn in monte_carlo_only_functions) {
      expect_false(
        grepl(paste0("\\b", mc_fn, "\\("), fn_src),
        info = paste(bu_fn, "should not call", mc_fn)
      )
    }
  }
})

test_that("structural boundary: fuse_observed_data_into_priors()/fuse_one_property_prior() are the ONLY monte-carlo.R functions bridging into bayesian-updating.R", {
  # The Bayesian-updating wiring (fuse_observed_data_into_priors(), and its
  # helper fuse_one_property_prior()) is the intentional, sole bridge between
  # the two files - every OTHER monte-carlo.R function should remain exactly
  # as decoupled from bayesian-updating.R as before this feature was added.
  bayesian_updating_exports <- c(
    "bayes_update_normal_normal", "normal_to_lognormal_params", "lognormal_to_normal_params",
    "fuse_beta", "fuse_gamma", "moments_to_beta", "moments_to_gamma", "beta_to_moments",
    "gamma_to_moments", "bayes_fuse", "bayesian_update", "fuse_bivariate_normal",
    "fuse_texture_group_from_triplets", "fuse_property"
  )
  other_monte_carlo_functions <- c(
    "generate_monte_carlo_realizations", "simulate_correlated_properties", "sim_component_compositions",
    "setup_distributions", "prepare_simulation_parameters", "configure_correlation_structure",
    "apply_simulation_constraints", "extract_property_parameters", "estimate_property_correlations",
    "prepare_simulation_data", "run_sequential_simulation", "run_parallel_simulation",
    "validate_monte_carlo_inputs", "validate_simulation_output", "get_monte_carlo_defaults",
    "normalize_monte_carlo_config", "get_constraint_rules", "get_sum_constraints",
    "as_lrh_triplet"
  )

  for (mc_fn in other_monte_carlo_functions) {
    fn_obj <- get(mc_fn, envir = asNamespace("soilSIM"))
    fn_src <- paste(deparse(body(fn_obj)), collapse = "\n")
    for (bu_fn in bayesian_updating_exports) {
      expect_false(
        grepl(paste0("\\b", bu_fn, "\\("), fn_src),
        info = paste(mc_fn, "should not call", bu_fn)
      )
    }
  }

  # Positive confirmation the bridge functions DO reference bayesian-updating.R,
  # so this test would actually fail (not just vacuously pass) if the bridge
  # were ever accidentally removed.
  bridge_src <- paste(
    deparse(body(get("fuse_observed_data_into_priors", envir = asNamespace("soilSIM")))),
    deparse(body(get("fuse_one_property_prior", envir = asNamespace("soilSIM")))),
    collapse = "\n"
  )
  expect_true(grepl("fuse_texture_group_from_triplets\\(", bridge_src))
  expect_true(grepl("bayes_update_normal_normal\\(", bridge_src))
})
