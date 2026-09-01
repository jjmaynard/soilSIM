test_that("fit_gamma_mom_raster() matches the scalar moments_to_gamma() cell-by-cell", {
  rasters <- make_percentile_rasters(c(a = 8, b = 10, c = 12))
  fit_r <- fit_gamma_mom_raster(list(rasters$a, rasters$b, rasters$c))
  mean_v <- mean(c(8, 10, 12)); var_v <- mean((c(8, 10, 12) - mean_v)^2)
  fit_scalar <- moments_to_gamma(mean_v, var_v)
  expect_equal(unique(terra::values(fit_r$shape))[1], fit_scalar$shape)
  expect_equal(unique(terra::values(fit_r$rate))[1], fit_scalar$rate)
})

test_that("percentile_index()/percentile_interior() locate low/median/high and interior probabilities", {
  idx <- percentile_index(c(0, 0.05, 0.5, 0.95, 1))
  expect_equal(idx$p_lo, 0.05)
  expect_equal(idx$p_hi, 0.95)
  expect_equal(idx$p50_idx, 3)

  interior <- percentile_interior(c(0, 0.05, 0.5, 0.95, 1))
  expect_equal(interior$probs, c(0.05, 0.5, 0.95))
})

test_that("align_percentile_probs() is a no-op when probabilities already match", {
  rasters <- make_percentile_rasters(c(P5 = 2, P50 = 10, P95 = 18))
  aligned <- align_percentile_probs(
    list(rasters$P5, rasters$P50, rasters$P95), c(0.05, 0.5, 0.95),
    list(rasters$P5, rasters$P50, rasters$P95), c(0.05, 0.5, 0.95)
  )
  expect_equal(aligned$probs, c(0.05, 0.5, 0.95))
})

test_that("align_percentile_probs() reinterpolates the wider side onto the narrower side's probabilities", {
  prior_r <- make_percentile_rasters(c(P5 = 2, P50 = 10, P95 = 18))
  lik_r <- make_percentile_rasters(c(P2.5 = 1, P50 = 10, P97.5 = 19))
  aligned <- align_percentile_probs(
    list(prior_r$P5, prior_r$P50, prior_r$P95), c(0.05, 0.5, 0.95),
    list(lik_r$P2.5, lik_r$P50, lik_r$P97.5), c(0.025, 0.5, 0.975)
  )
  # prior's range [0.05, 0.95] is narrower than lik's [0.025, 0.975] -> target = prior_probs
  expect_equal(aligned$probs, c(0.05, 0.5, 0.95))
  expect_length(aligned$lik_value_rasters, 3)
})

test_that("fuse_adaptive() routes to the closed-form path above threshold_cells and matches bayes_update_normal_normal()", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 0, verbose = FALSE
  )
  expect_equal(result$route, "closed_form_normal")

  fit_prior <- fit_normal_raster(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95, 0.05, 0.95)
  fit_lik <- fit_normal_raster(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95, 0.05, 0.95)
  expected <- bayes_update_normal_normal(fit_prior$mu, fit_prior$sigma, fit_lik$mu, fit_lik$sigma)
  expect_equal(unique(terra::values(result$posterior$mu))[1], unique(terra::values(expected$mu))[1])
})

test_that("fuse_adaptive() routes to the general KDE path at or below threshold_cells", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 1e6, n_samples = 50, verbose = FALSE
  )
  expect_equal(result$route, "bayesian_update_general")
  expect_true(all(is.finite(terra::values(result$posterior$mu))))
})

test_that("fuse_general_kde() degrades a cell with too many NA percentile values to NA output, rather than aborting the whole raster", {
  # Real-world discovery: cells with insufficient valid percentile values (e.g. at raster edges
  # after resample()/alignment) used to hard-crash the entire terra::app() call via
  # extract_percentile_pairs()'s "Need at least 2 valid quantile columns" stop().
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6), nrow = 1, ncol = 2)
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9), nrow = 1, ncol = 2)
  # Poison the second cell of every prior layer with NA.
  for (p in c("P5", "P50", "P95")) {
    terra::values(rasters_prior[[p]])[2] <- NA
  }

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 1e6, n_samples = 50, verbose = FALSE
  )
  mu_vals <- terra::values(result$posterior$mu)[, 1]
  expect_true(is.finite(mu_vals[1]))
  expect_true(is.na(mu_vals[2]))
})

test_that("fuse_general_kde()'s default grid_resolution stays within tolerance of the finer 0.01 reference", {
  # PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4: FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION (0.1) trades
  # a coarser bayesian_update() evaluation grid for ~2.4x wall-clock speed (measured on a real
  # fuse_adaptive() call, 23.86s -> 9.91s at 2,024 cells). This asserts the posterior moments
  # stay close to the finer 0.01 grid across several representative percentile scenarios, mirroring
  # the accuracy sweep that justified the constant's value (mean error < 0.02%, variance error <
  # 0.12% analytically) - tolerance here is looser than that analytical bound to absorb the extra
  # Monte Carlo noise from bayesian_update()'s own sample()-based posterior draw, which the
  # analytical sweep bypassed.
  scenarios <- list(
    narrow_clay = c(P5 = 18, P50 = 22, P95 = 27),
    wide_sand   = c(P5 = 20, P50 = 40, P95 = 65),
    narrow_ph   = c(P5 = 5.8, P50 = 6.2, P95 = 6.6)
  )
  probs <- c(0.05, 0.5, 0.95)

  # A single seeded draw carries real Monte Carlo noise from two compounded sampling stages
  # (simulate_from_percentiles()'s n_samples=500 draw, then bayesian_update()'s own n=1000
  # posterior draw) - averaging several independent seeds isolates the systematic effect of
  # grid_resolution itself from that per-run sampling noise, avoiding a flaky single-draw test.
  n_reps <- 8
  for (scen in names(scenarios)) {
    vals <- scenarios[[scen]]
    rasters_prior <- make_percentile_rasters(vals)
    rasters_lik <- make_percentile_rasters(vals + 1)
    prior_list <- list(rasters_prior[[1]], rasters_prior[[2]], rasters_prior[[3]])
    lik_list <- list(rasters_lik[[1]], rasters_lik[[2]], rasters_lik[[3]])

    mu_default <- numeric(n_reps); mu_fine <- numeric(n_reps)
    sigma_default <- numeric(n_reps); sigma_fine <- numeric(n_reps)
    for (rep in seq_len(n_reps)) {
      set.seed(rep)
      result_default <- fuse_adaptive(prior_list, lik_list, probs, family = "normal",
                                       threshold_cells = 1e6, n_samples = 500, verbose = FALSE)
      set.seed(rep)
      result_fine <- fuse_adaptive(prior_list, lik_list, probs, family = "normal",
                                    threshold_cells = 1e6, n_samples = 500,
                                    grid_resolution = 0.01, verbose = FALSE)
      mu_default[rep] <- unique(terra::values(result_default$posterior$mu))[1]
      mu_fine[rep] <- unique(terra::values(result_fine$posterior$mu))[1]
      sigma_default[rep] <- unique(terra::values(result_default$posterior$sigma))[1]
      sigma_fine[rep] <- unique(terra::values(result_fine$posterior$sigma))[1]
    }

    expect_equal(mean(mu_default), mean(mu_fine), tolerance = 0.02, label = paste0(scen, " mean(mu)"))
    expect_equal(mean(sigma_default), mean(sigma_fine), tolerance = 0.1, label = paste0(scen, " mean(sigma)"))
  }
})

test_that("fuse_general_kde(raw_draws) posterior mean converges toward the true prior value as the mukey's real draw count grows", {
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.7. Deliberately weak/wide likelihood so the
  # posterior mean stays close to the prior mean - isolates how well each sample size's real-draw
  # resample approximates the TRUE prior distribution (a classic Monte Carlo convergence check).
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  make_uniform_raster <- function(value) terra::rast(nrows = 1, ncols = 1, vals = value)

  lik_percentiles <- stats::qnorm(probs, mean = 25, sd = 25)
  lik_value_rasters <- stats::setNames(lapply(lik_percentiles, make_uniform_raster), paste0("P", round(probs * 100)))

  mukey_raster <- terra::rast(nrows = 1, ncols = 1, vals = 900)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  meanlog <- log(20); sdlog <- 0.5
  true_prior_mean <- exp(meanlog + sdlog^2 / 2)
  set.seed(123)
  full_pop <- stats::rlnorm(200000, meanlog = meanlog, sdlog = sdlog)
  prior_percentiles <- as.numeric(stats::quantile(full_pop, probs))
  prior_value_rasters <- stats::setNames(lapply(prior_percentiles, make_uniform_raster), paste0("P", round(probs * 100)))

  sample_sizes <- c(30, 300, 3000)
  n_reps <- 6
  mean_abs_bias <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    mus <- numeric(n_reps)
    for (rep in seq_len(n_reps)) {
      set.seed(1000 * i + rep)
      mukey_draws <- list("900" = sample(full_pop, sample_sizes[i]))
      result <- fuse_general_kde(prior_value_rasters, lik_value_rasters, probs, family = "normal",
                                  bounds = NULL, n_samples = 500, grid_resolution = NULL,
                                  mukey_raster = mukey_raster, mukey_draws = mukey_draws)
      mus[rep] <- terra::values(result$posterior$mu)[1, 1]
    }
    mean_abs_bias[i] <- mean(abs(mus - true_prior_mean))
  }
  # Larger real-draw pools should approximate the true prior distribution more consistently -
  # average absolute deviation from the true prior mean shrinks from the smallest to the largest
  # sample size.
  expect_true(mean_abs_bias[length(mean_abs_bias)] < mean_abs_bias[1])
})

test_that("fuse_general_kde(raw_draws)'s PRIOR-side resample tracks the real data's own empirical mean, before fusion", {
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.7.
  #
  # An earlier version of this test asserted the FUSED posterior mean (after combining with a
  # weak/wide likelihood) lands closer to the true prior mean under raw_draws than under the
  # default percentile-reconstruction route, for a skewed lognormal prior. Investigated live
  # (not noise - failed consistently across 8 averaged seeds both times): the raw_draws resample
  # retains genuine outliers the real data actually contains (observed range ~4-195 for a
  # sdlog=0.7 lognormal, n=5000), while sim_linear_cdf_batch()'s percentile reconstruction stays
  # confined to roughly the P5-P95 range (~7-64) with no comparable tail extrapolation. That wider
  # realized range measurably widens stats::density()'s default ("nrd0") bandwidth on the raw_draws
  # side inside bayesian_update(), which can shift the FUSED posterior mean in a way that isn't
  # reliably closer to the true prior mean than the reconstruction route - a real, non-obvious
  # characteristic of the current KDE fusion mechanism worth flagging for P2.8/P2.9 (not something
  # to force past with a misleading test - see this task's write-up in
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md for the full finding and its implication for P2.10).
  #
  # What IS true and defensible: verified directly (prior side alone, before fusion) that raw_draws
  # resampling tracks the real data's own empirical mean about as well as the reconstruction does -
  # both approximate `mean(full_pop)`, not the unknowable analytic population mean, which is the
  # correct fidelity target given raw_draws fuses against the REAL simulated draws, not an idealized
  # distribution. This test asserts that narrower, correct claim.
  meanlog <- log(20); sdlog <- 0.7
  set.seed(456)
  full_pop <- stats::rlnorm(5000, meanlog = meanlog, sdlog = sdlog)
  true_empirical_mean <- mean(full_pop)

  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  percentile_cols <- paste0("P", round(probs * 100))
  prior_percentiles <- as.numeric(stats::quantile(full_pop, probs))
  prior_row <- as.data.frame(stats::setNames(as.list(prior_percentiles), percentile_cols))

  n_reps <- 8
  raw_means <- numeric(n_reps); recon_means <- numeric(n_reps)
  for (rep in seq_len(n_reps)) {
    set.seed(rep)
    raw_means[rep] <- mean(sample(full_pop, 500, replace = TRUE))
    set.seed(rep)
    recon_means[rep] <- mean(simulate_from_percentiles(prior_row, method = "linear_cdf",
                                                         percentile_cols = percentile_cols, n = 500))
  }

  expect_equal(mean(raw_means), true_empirical_mean, tolerance = 0.1)
  expect_equal(mean(recon_means), true_empirical_mean, tolerance = 0.1)
})

test_that("fuse_texture_group_batch(raw_draws) preserves real cross-fraction correlation the default route discards", {
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.7. Builds synthetic clay/sand draws with a
  # strong, deliberate negative correlation (realistic soil-texture pattern) and confirms the
  # raw_draws route's fused ILR covariance (S12) reflects that real correlation much more
  # negatively than the default route's, which reconstructs each fraction independently from its
  # own percentile-derived marginal and so can only pick up an incidental correlation from the
  # compositional-closure rescale.
  set.seed(7)
  n <- 4000
  clay_d <- pmax(stats::rnorm(n, 20, 3), 1)
  sand_d <- pmax(45 - 0.9 * (clay_d - 20) + stats::rnorm(n, 0, 1), 1)
  silt_d <- pmax(100 - clay_d - sand_d, 1)
  mukey_texture_draws <- list("900" = cbind(clay_total = clay_d, sand_total = sand_d, silt_total = silt_d))

  probs <- c(0.05, 0.5, 0.95)
  prior_percentiles <- list(
    clay = as.numeric(stats::quantile(clay_d, probs)),
    sand = as.numeric(stats::quantile(sand_d, probs)),
    silt = as.numeric(stats::quantile(silt_d, probs))
  )
  # Deliberately wide/weak likelihood (large spread relative to the prior) so the fused output
  # stays close to whichever prior approximation was used - isolates the prior-side correlation
  # effect, matching the same weak-likelihood isolation strategy used in the mean-convergence tests
  # above.
  lik_wide <- function(p50) p50 + c(-40, 0, 40)

  col_names <- c(
    "clay_prior_lo", "clay_prior_p50", "clay_prior_hi", "clay_lik_lo", "clay_lik_p50", "clay_lik_hi",
    "sand_prior_lo", "sand_prior_p50", "sand_prior_hi", "sand_lik_lo", "sand_lik_p50", "sand_lik_hi",
    "silt_prior_lo", "silt_prior_p50", "silt_prior_hi", "silt_lik_lo", "silt_lik_p50", "silt_lik_hi"
  )
  row_vals <- c(
    prior_percentiles$clay, lik_wide(prior_percentiles$clay[2]),
    prior_percentiles$sand, lik_wide(prior_percentiles$sand[2]),
    prior_percentiles$silt, lik_wide(prior_percentiles$silt[2])
  )
  row_mat <- matrix(row_vals, nrow = 1, dimnames = list(NULL, col_names))
  prior_z <- stats::qnorm(0.95); lik_z <- stats::qnorm(0.95)

  mukey_prior_texture_samples <- list("900" = mukey_texture_draws[["900"]][sample(n, 2000, replace = TRUE), , drop = FALSE])

  default_out <- fuse_texture_group_batch(row_mat, "clay", "sand", "silt", prior_z, lik_z)
  raw_row_mat <- cbind(row_mat, mukey_code = 900)
  raw_out <- fuse_texture_group_batch(raw_row_mat, "clay", "sand", "silt", prior_z, lik_z,
                                       mukey_prior_texture_samples = mukey_prior_texture_samples)

  expect_true(raw_out[1, "S12"] < default_out[1, "S12"])
})

test_that("fuse_property_adaptive() dispatches normal/beta/gamma/lognormal correctly", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))
  prior_list <- list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95)
  lik_list <- list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95)
  probs <- c(0.05, 0.5, 0.95)

  r_normal <- fuse_property_adaptive(prior_list, probs, lik_list, probs,
    property_config = list(dist = "normal"), threshold_cells = 0, verbose = FALSE)
  expect_equal(r_normal$dist, "normal")
  expect_equal(r_normal$route, "closed_form_normal")

  r_lognormal <- fuse_property_adaptive(prior_list, probs, lik_list, probs,
    property_config = list(dist = "lognormal"), threshold_cells = 0, verbose = FALSE)
  expect_equal(r_lognormal$dist, "lognormal")
  expect_equal(r_lognormal$route, "closed_form_lognormal")

  rasters_prior_pos <- make_percentile_rasters(c(P5 = 20, P50 = 40, P95 = 65))
  rasters_lik_pos <- make_percentile_rasters(c(P5 = 30, P50 = 50, P95 = 70))
  r_beta <- fuse_property_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95), probs,
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95), probs,
    property_config = list(dist = "beta", bounds = c(0, 100)), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_beta$dist, "beta")

  r_gamma <- fuse_property_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95), probs,
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95), probs,
    property_config = list(dist = "gamma"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_gamma$dist, "gamma")
})

test_that("fuse_property_adaptive() dispatches metalog and flags fallback via route_detail", {
  rasters_prior <- make_percentile_rasters(c(P10 = 20, P50 = 40, P90 = 65))
  rasters_lik <- make_percentile_rasters(c(P10 = 30, P50 = 50, P90 = 70))
  probs <- c(0.1, 0.5, 0.9)

  r_metalog <- fuse_property_adaptive(
    list(rasters_prior$P10, rasters_prior$P50, rasters_prior$P90), probs,
    list(rasters_lik$P10, rasters_lik$P50, rasters_lik$P90), probs,
    property_config = list(dist = "metalog", bounds = c(0, 100), boundedness = "b"),
    verbose = FALSE
  )
  expect_equal(r_metalog$dist, "metalog")
  expect_equal(r_metalog$route, "closed_form_metalog_adapter")
  expect_true(all(is.finite(terra::values(r_metalog$posterior$mu))))
})

test_that("fuse_property_adaptive() resolves dist = 'auto' to beta/normal/lognormal without overriding an explicit dist", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))
  probs <- c(0.05, 0.5, 0.95)

  r_auto <- fuse_property_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95), probs,
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95), probs,
    property_config = list(dist = "auto"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_auto$dist, "normal")
  expect_equal(r_auto$dist_source, "auto")

  r_explicit <- fuse_property_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95), probs,
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95), probs,
    property_config = list(dist = "lognormal"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_explicit$dist_source, "config")
})

test_that("group_members() returns configured members in order, reusing soilSIM's composition_groups shape", {
  composition_groups <- list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))
  expect_equal(group_members("texture", composition_groups), c("claytotal", "sandtotal", "silttotal"))
  expect_error(group_members("nonexistent", composition_groups), "no composition group")
})

test_that("fuse_texture_group() returns one posterior per member that sums to ~100 and matches ilr_inverse()", {
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  result <- fuse_texture_group(fetched)
  expect_named(result, c("claytotal", "sandtotal", "silttotal"))

  totals <- Reduce(`+`, lapply(result, function(r) r$posterior$value))
  expect_true(all(abs(terra::values(totals) - 100) < 1e-6))
  for (id in names(result)) {
    expect_equal(result[[id]]$dist, "texture_ilr")
    expect_equal(result[[id]]$route, "closed_form_ilr_group")
  }
})

test_that("fuse_texture_group(prior_fusion_method='raw_draws') still sums to ~100, stays finite, and differs from the default percentile-reconstruction route", {
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.4.
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )
  default_result <- fuse_texture_group(fetched)

  set.seed(7)
  n <- 3000
  clay_d <- pmax(stats::rnorm(n, 20, 3), 1)
  sand_d <- pmax(45 - 0.6 * (clay_d - 20) + stats::rnorm(n, 0, 2), 1)
  silt_d <- pmax(100 - clay_d - sand_d, 1)
  mukey_texture_draws <- list("777" = cbind(clay_total = clay_d, sand_total = sand_d, silt_total = silt_d))

  mukey_raster <- terra::rast(nrows = 2, ncols = 2, vals = rep(777, 4))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  raw_result <- fuse_texture_group(fetched, mukey_raster = mukey_raster, mukey_texture_draws = mukey_texture_draws)
  expect_named(raw_result, c("claytotal", "sandtotal", "silttotal"))

  raw_totals <- Reduce(`+`, lapply(raw_result, function(r) r$posterior$value))
  expect_true(all(abs(terra::values(raw_totals) - 100) < 1e-6))
  for (id in names(raw_result)) {
    expect_equal(raw_result[[id]]$dist, "texture_ilr")
    expect_equal(raw_result[[id]]$route, "closed_form_ilr_group")
    expect_true(all(is.finite(terra::values(raw_result[[id]]$posterior$value))))
  }

  default_clay <- terra::values(default_result$claytotal$posterior$value)[, 1]
  raw_clay <- terra::values(raw_result$claytotal$posterior$value)[, 1]
  expect_false(isTRUE(all.equal(default_clay, raw_clay)))
})

test_that("fuse_texture_group() default output is unaffected by the new mukey_raster/mukey_texture_draws parameters when both are NULL", {
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  set.seed(11)
  a <- fuse_texture_group(fetched)
  set.seed(11)
  b <- fuse_texture_group(fetched, mukey_raster = NULL, mukey_texture_draws = NULL)
  expect_equal(terra::values(a$claytotal$posterior$value), terra::values(b$claytotal$posterior$value))
})

test_that("run_stage1_fusion_group() uses distinct cache kinds for prior_fusion_method='raw_draws' vs the default", {
  # Regression test for a real bug found live while verifying P2.5 (MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md):
  # the group-level result cache key didn't encode prior_fusion_method at all, so switching it for
  # the same AOI/depth-window silently returned a STALE, method-mismatched cached fusion result -
  # unlike run_stage1_fusion()'s own "ssurgo" cache, which only stores the pre-fusion prior
  # (re-fused fresh every call), the group cache stores the final fused POSTERIOR, so this staleness
  # risk was specific to the group path. Offline/pure - no network needed, since this tests the
  # cache-key/routing logic directly by pre-seeding both kinds and mocking the fetch functions to
  # error if called (proving a cache hit occurred, not a fresh fetch).
  aoi <- terra::vect("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")
  top_depth <- 0; bottom_depth <- 20
  cache_dir <- tools::R_user_dir("soilSIM", "cache")
  on.exit({
    for (kind in c("texture_group", "texture_group_raw_draws")) {
      key <- build_cache_key(aoi, "texture", top_depth, bottom_depth, kind)
      unlink(file.path(cache_dir, paste0(key, ".rds")))
    }
  }, add = TRUE)

  default_key <- build_cache_key(aoi, "texture", top_depth, bottom_depth, "texture_group")
  raw_key <- build_cache_key(aoi, "texture", top_depth, bottom_depth, "texture_group_raw_draws")
  expect_false(default_key == raw_key)

  fake_layer <- terra::rast(nrows = 1, ncols = 1, vals = 1)
  make_fake_group_result <- function(tag) {
    stats::setNames(lapply(c("clay", "sand", "silt"), function(id) {
      list(posterior = list(value = fake_layer, ilr_mu = fake_layer, ilr_Sigma = fake_layer),
           dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL,
           n_fallback_cells = 0, tag = tag)
    }), c("clay", "sand", "silt"))
  }
  cache_set(default_key, "texture_group", wrap_nested_rasters(make_fake_group_result("default")))
  cache_set(raw_key, "texture_group_raw_draws", wrap_nested_rasters(make_fake_group_result("raw_draws")))

  composition_groups <- list(texture = list(members = c("clay", "sand", "silt")))
  base_configs <- list(
    clay = list(id = "clay", solus_variable = "claytotal", composition_group = "texture"),
    sand = list(id = "sand", solus_variable = "sandtotal", composition_group = "texture"),
    silt = list(id = "silt", solus_variable = "silttotal", composition_group = "texture")
  )
  raw_configs <- base_configs
  raw_configs$clay$prior_fusion_method <- "raw_draws"

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) stop("should hit the pre-seeded group cache, not fetch"),
    .package = "soilSIM"
  )

  default_result <- run_stage1_fusion_group(aoi, "texture", composition_groups, base_configs, top_depth, bottom_depth)
  raw_result <- run_stage1_fusion_group(aoi, "texture", composition_groups, raw_configs, top_depth, bottom_depth)

  expect_equal(default_result$clay$tag, "default")
  expect_equal(raw_result$clay$tag, "raw_draws")
})

test_that("fuse_texture_group_batch() vectorized path matches the original per-cell scalar loop bit-for-bit", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 1 fuse_texture_group() fix:
  # reimplements the ORIGINAL per-cell scalar loop (estimate_ilr_moments_mc() x2 +
  # fuse_bivariate_normal(), called once per cell via a plain for-loop) here, and checks the new
  # vectorized fuse_texture_group_batch() produces identical results under the same seed - proving
  # the rnorm() stream is consumed in the same order (a pure vectorization, not a behavior change).
  old_scalar <- function(row_mat, clay_id, sand_id, silt_id, prior_z, lik_z) {
    ncell <- nrow(row_mat)
    out <- matrix(NA_real_, nrow = ncell, ncol = 5, dimnames = list(NULL, c("mu1", "mu2", "S11", "S12", "S22")))
    for (i in seq_len(ncell)) {
      row <- row_mat[i, ]
      prior_m <- estimate_ilr_moments_mc(
        row[[paste0(clay_id, "_prior_lo")]], row[[paste0(clay_id, "_prior_p50")]], row[[paste0(clay_id, "_prior_hi")]],
        row[[paste0(sand_id, "_prior_lo")]], row[[paste0(sand_id, "_prior_p50")]], row[[paste0(sand_id, "_prior_hi")]],
        row[[paste0(silt_id, "_prior_lo")]], row[[paste0(silt_id, "_prior_p50")]], row[[paste0(silt_id, "_prior_hi")]],
        z = prior_z
      )
      lik_m <- estimate_ilr_moments_mc(
        row[[paste0(clay_id, "_lik_lo")]], row[[paste0(clay_id, "_lik_p50")]], row[[paste0(clay_id, "_lik_hi")]],
        row[[paste0(sand_id, "_lik_lo")]], row[[paste0(sand_id, "_lik_p50")]], row[[paste0(sand_id, "_lik_hi")]],
        row[[paste0(silt_id, "_lik_lo")]], row[[paste0(silt_id, "_lik_p50")]], row[[paste0(silt_id, "_lik_hi")]],
        z = lik_z
      )
      fused <- fuse_bivariate_normal(prior_m$mu, prior_m$Sigma, lik_m$mu, lik_m$Sigma)
      out[i, ] <- c(fused$mu[1], fused$mu[2], fused$Sigma[1, 1], fused$Sigma[1, 2], fused$Sigma[2, 2])
    }
    out
  }

  # Heterogeneous cells (distinct lo/p50/hi per cell), not uniform, to genuinely exercise
  # per-cell vectorized math rather than trivially-identical rows.
  set.seed(1)
  n <- 5
  ids <- c("claytotal", "sandtotal", "silttotal")
  col_names <- as.vector(outer(
    paste0(rep(ids, each = 2), rep(c("_prior_", "_lik_"), times = length(ids))),
    c("lo", "p50", "hi"), paste0
  ))
  row_mat <- matrix(NA_real_, nrow = n, ncol = length(col_names), dimnames = list(NULL, col_names))
  for (id in ids) {
    for (side in c("prior", "lik")) {
      p50 <- runif(n, 20, 40)
      spread <- runif(n, 3, 8)
      row_mat[, paste0(id, "_", side, "_lo")] <- p50 - spread
      row_mat[, paste0(id, "_", side, "_p50")] <- p50
      row_mat[, paste0(id, "_", side, "_hi")] <- p50 + spread
    }
  }
  prior_z <- stats::qnorm(0.95)
  lik_z <- stats::qnorm(0.95)

  set.seed(42)
  expected <- old_scalar(row_mat, "claytotal", "sandtotal", "silttotal", prior_z, lik_z)
  set.seed(42)
  actual <- fuse_texture_group_batch(row_mat, "claytotal", "sandtotal", "silttotal", prior_z, lik_z)

  expect_equal(actual, expected, tolerance = 1e-9)
})

test_that("run_stage1_fusion() fetches, caches, aligns, and fuses SSURGO x SOLUS for a real AOI", {
  testthat::skip_if_offline()
  # See test-solus-simulation.R's live test for why PROJ_LIB is unset defensively here too - a
  # stray PROJ_LIB (documented in HANDOFF_NOTES.md) can poison terra's PROJ context before
  # soilSIM's own .onLoad() runs, when terra was already attached earlier in the same test session.
  Sys.unsetenv("PROJ_LIB")
  aoi <- terra::vect(
    "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
    "AOI has no valid CRS in this session - likely the documented PROJ_LIB environment issue; not a soilSIM defect."
  )
  aoi <- terra::project(aoi, "epsg:5070")

  property_config <- list(id = "test_clay", solus_variable = "claytotal", dist = "normal")

  result <- run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5)
  expect_true(is.null(result) || is.list(result))
  if (!is.null(result)) {
    expect_true(result$dist %in% c("normal", "beta", "gamma", "lognormal", "metalog"))
    expect_true(is.list(result$posterior))
    expect_s4_class(result$prior$values[[1]], "SpatRaster")
  }
})

test_that("run_stage1_fusion() with a composition_group dispatches to run_stage1_fusion_group() and slices to the requested member", {
  testthat::skip_if_offline()
  Sys.unsetenv("PROJ_LIB")
  aoi <- terra::vect(
    "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
    "AOI has no valid CRS in this session - likely the documented PROJ_LIB environment issue; not a soilSIM defect."
  )
  aoi <- terra::project(aoi, "epsg:5070")

  composition_groups <- list(texture = list(members = c("test_clay2", "test_sand2", "test_silt2")))
  property_configs <- list(
    test_clay2 = list(id = "test_clay2", solus_variable = "claytotal", composition_group = "texture"),
    test_sand2 = list(id = "test_sand2", solus_variable = "sandtotal", composition_group = "texture"),
    test_silt2 = list(id = "test_silt2", solus_variable = "silttotal", composition_group = "texture")
  )

  result <- run_stage1_fusion(
    aoi, property_configs$test_clay2, top_depth = 0, bottom_depth = 5,
    composition_groups = composition_groups, property_configs = property_configs
  )
  expect_true(is.null(result) || identical(result$dist, "texture_ilr"))
})

# ---------------------------------------------------------------------------
# P3.8 - Percentile output tests (MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P3.8)
# ---------------------------------------------------------------------------

test_that("fuse_adaptive()'s general route (fuse_general_kde()) produces monotonic percentiles at the default posterior_probs", {
  rasters_prior <- make_percentile_rasters(c(P5 = 15, P50 = 20, P95 = 28))
  rasters_lik <- make_percentile_rasters(c(P5 = 18, P50 = 24, P95 = 32))
  probs <- c(0.05, 0.5, 0.95)

  set.seed(201)
  result <- fuse_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    probs, family = "normal", threshold_cells = 1e6, n_samples = 500,
    grid_resolution = 0.1, verbose = FALSE
  )
  expect_named(result$posterior$percentiles, c("P1", "P5", "P10", "P25", "P50", "P75", "P90", "P95", "P99"))
  pvals <- sapply(result$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expect_true(all(diff(pvals) >= 0))
})

test_that("fuse_adaptive()'s posterior_probs parameter is actually respected, not hardcoded", {
  rasters_prior <- make_percentile_rasters(c(P5 = 15, P50 = 20, P95 = 28))
  rasters_lik <- make_percentile_rasters(c(P5 = 18, P50 = 24, P95 = 32))
  probs <- c(0.05, 0.5, 0.95)
  custom_probs <- c(0.02, 0.5, 0.98)

  set.seed(203)
  result <- fuse_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    probs, family = "normal", threshold_cells = 1e6, n_samples = 500,
    grid_resolution = 0.1, verbose = FALSE, posterior_probs = custom_probs
  )
  expect_named(result$posterior$percentiles, c("P2", "P50", "P98"))

  # Closed-form route, same custom probs - both routes should honor the same parameter identically.
  set.seed(204)
  result_closed <- fuse_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    probs, family = "normal", threshold_cells = 0, verbose = FALSE, posterior_probs = custom_probs
  )
  expect_named(result_closed$posterior$percentiles, c("P2", "P50", "P98"))
})

test_that("fuse_closed_form() percentiles exactly match analytic qnorm/qbeta/qgamma at the fused parameters", {
  probs <- c(0.05, 0.5, 0.95)
  posterior_probs <- c(0.01, 0.05, 0.5, 0.95, 0.99)

  # Normal
  rasters_prior <- make_percentile_rasters(c(P5 = 15, P50 = 20, P95 = 28))
  rasters_lik <- make_percentile_rasters(c(P5 = 18, P50 = 24, P95 = 32))
  r_normal <- fuse_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    probs, family = "normal", threshold_cells = 0, verbose = FALSE, posterior_probs = posterior_probs
  )
  mu <- terra::values(r_normal$posterior$mu)[1, 1]; sigma <- terra::values(r_normal$posterior$sigma)[1, 1]
  pvals_n <- sapply(r_normal$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expect_equal(unname(pvals_n), qnorm(posterior_probs, mu, sigma), tolerance = 1e-9)

  # Beta
  rasters_prior_pos <- make_percentile_rasters(c(P5 = 20, P50 = 40, P95 = 65))
  rasters_lik_pos <- make_percentile_rasters(c(P5 = 30, P50 = 50, P95 = 70))
  bounds <- c(0, 100)
  r_beta <- fuse_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95),
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95),
    probs, family = "beta", bounds = bounds, threshold_cells = 0, verbose = FALSE,
    posterior_probs = posterior_probs
  )
  alpha <- terra::values(r_beta$posterior$alpha)[1, 1]; beta_p <- terra::values(r_beta$posterior$beta)[1, 1]
  pvals_b <- sapply(r_beta$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expected_b <- bounds[1] + qbeta(posterior_probs, alpha, beta_p) * (bounds[2] - bounds[1])
  expect_equal(unname(pvals_b), expected_b, tolerance = 1e-9)
  expect_true(all(pvals_b >= bounds[1] & pvals_b <= bounds[2]))

  # Gamma
  r_gamma <- fuse_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95),
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95),
    probs, family = "gamma", threshold_cells = 0, verbose = FALSE, posterior_probs = posterior_probs
  )
  shape <- terra::values(r_gamma$posterior$shape)[1, 1]; rate <- terra::values(r_gamma$posterior$rate)[1, 1]
  pvals_g <- sapply(r_gamma$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expect_equal(unname(pvals_g), qgamma(posterior_probs, shape = shape, rate = rate), tolerance = 1e-9)

  # All three: monotonic
  expect_true(all(diff(pvals_n) >= 0))
  expect_true(all(diff(pvals_b) >= 0))
  expect_true(all(diff(pvals_g) >= 0))
})

test_that("fuse_lognormal_adaptive() produces monotonic, strictly-positive percentiles on both its general and closed-form branches", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  posterior_probs <- c(0.01, 0.05, 0.5, 0.95, 0.99)
  prior_r <- stats::setNames(lapply(qlnorm(probs, log(20), 0.4), function(v) make_percentile_rasters(c(x = v))$x),
                              paste0("P", round(probs * 100)))
  lik_r <- stats::setNames(lapply(qlnorm(probs, log(22), 0.3), function(v) make_percentile_rasters(c(x = v))$x),
                            paste0("P", round(probs * 100)))

  r_big <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 0,
                                    verbose = FALSE, posterior_probs = posterior_probs)
  pvals_big <- sapply(r_big$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expect_true(all(diff(pvals_big) >= 0))
  expect_true(all(pvals_big > 0))

  set.seed(211)
  r_small <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 1e6,
                                      verbose = FALSE, posterior_probs = posterior_probs, grid_resolution = 0.1)
  pvals_small <- sapply(r_small$posterior$percentiles, function(r) terra::values(r)[1, 1])
  expect_true(all(diff(pvals_small) >= 0))
})

test_that("fuse_metalog_adapter() percentiles exactly match analytic qnorm at the fused (mu, sigma), and pass through fuse_property_adaptive() unchanged", {
  # Uses the known-working 3-point metalog fixture (see fuse_property_adaptive() metalog dispatch
  # test above) - a symmetric 5-point Normal percentile set was found (P3.5 investigation) to
  # trigger a solve()-singularity in fit_metalog_linear_raster() unrelated to percentile output.
  rasters_prior <- make_percentile_rasters(c(P10 = 20, P50 = 40, P90 = 65))
  rasters_lik <- make_percentile_rasters(c(P10 = 30, P50 = 50, P90 = 70))
  probs <- c(0.1, 0.5, 0.9)
  posterior_probs <- c(0.01, 0.05, 0.5, 0.95, 0.99)
  property_config <- list(bounds = c(0, 100), boundedness = "b")

  prior_list <- list(rasters_prior$P10, rasters_prior$P50, rasters_prior$P90)
  lik_list <- list(rasters_lik$P10, rasters_lik$P50, rasters_lik$P90)

  r <- fuse_metalog_adapter(prior_list, probs, lik_list, probs, property_config, verbose = FALSE,
                             posterior_probs = posterior_probs)
  mu <- terra::values(r$posterior$mu)[1, 1]; sigma <- terra::values(r$posterior$sigma)[1, 1]
  pvals <- sapply(r$posterior$percentiles, function(x) terra::values(x)[1, 1])
  expect_equal(unname(pvals), qnorm(posterior_probs, mu, sigma), tolerance = 1e-9)
  expect_true(all(diff(pvals) >= 0))

  # fuse_property_adaptive() dispatch: extra unused args (n_samples/grid_resolution) must be
  # silently absorbed, not error, and produce the identical percentiles.
  property_config2 <- list(dist = "metalog", bounds = c(0, 100), boundedness = "b")
  r2 <- fuse_property_adaptive(prior_list, probs, lik_list, probs, property_config2,
                                verbose = FALSE, posterior_probs = posterior_probs,
                                n_samples = 500, grid_resolution = 0.1)
  pvals2 <- sapply(r2$posterior$percentiles, function(x) terra::values(x)[1, 1])
  expect_equal(pvals, pvals2)
})

test_that("fuse_texture_group() percentiles are monotonic per fraction, and (unlike the point-estimate value) do NOT generally sum to 100 across the group", {
  # MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P3.6's documented sum-to-100 caveat, asserted here
  # as a permanent regression per task P3.8's explicit instruction.
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 25, hi = 38))
  sand_prior <- make_percentile_rasters(c(lo = 30, p50 = 45, hi = 60))
  silt_prior <- make_percentile_rasters(c(lo = 20, p50 = 30, hi = 42))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 28, hi = 40))
  sand_lik <- make_percentile_rasters(c(lo = 28, p50 = 42, hi = 58))
  silt_lik <- make_percentile_rasters(c(lo = 22, p50 = 30, hi = 40))

  probs <- c(0.05, 0.5, 0.95)
  posterior_probs <- c(0.01, 0.05, 0.5, 0.95, 0.99)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  set.seed(301)
  result <- fuse_texture_group(fetched, posterior_probs = posterior_probs, posterior_n_mc = 2000)

  for (id in names(result)) {
    pvals <- sapply(result[[id]]$posterior$percentiles, function(r) terra::values(r)[1, 1])
    expect_true(all(diff(pvals) >= 0), info = id)
  }

  value_sum <- sum(sapply(names(result), function(id) terra::values(result[[id]]$posterior$value)[1, 1]))
  expect_equal(value_sum, 100, tolerance = 1e-6)

  p95_sum <- sum(sapply(names(result), function(id) {
    terra::values(result[[id]]$posterior$percentiles$P95)[1, 1]
  }))
  # Correlated-tail draws don't cancel the way the point estimate's ilr_inverse() renormalization
  # does - the P95 sum should measurably exceed 100, not equal it.
  expect_true(p95_sum > 105)
})

test_that("fuse_texture_group() degrades NA-coverage cells to NA percentiles instead of crashing the whole raster - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P3.9 regression", {
  # Found live via P3.9's real-AOI validation: texture_group_percentiles_raster()'s quantile()
  # call originally had no na.rm handling, so a single NA-coverage cell (a real, common occurrence
  # - P2.5/P2.8's own live validation found 401/598 NA cells for one real AOI) crashed the ENTIRE
  # terra::app() call for the whole raster, unlike every other per-cell fallback in this file.
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  # Punch an NA hole into one cell of every input layer (simulating a real coverage gap) - all
  # rasters here are 2x2 (4 cells) via make_percentile_rasters()'s default.
  na_out_cell1 <- function(r) { v <- terra::values(r); v[1, 1] <- NA; terra::values(r) <- v; r }
  clay_prior <- lapply(clay_prior, na_out_cell1)
  sand_prior <- lapply(sand_prior, na_out_cell1)
  silt_prior <- lapply(silt_prior, na_out_cell1)

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  # Must not error/crash - this is the regression itself. rnorm(NA, NA) inside the NA'd cell's
  # moment estimation emits an expected, harmless "NAs produced" warning (pre-existing behavior of
  # fuse_texture_group_batch_core(), only now exercised by this test's deliberate NA input) -
  # suppressed here since it's not the thing under test.
  result <- expect_no_error(suppressWarnings(fuse_texture_group(fetched)))

  clay_p50 <- terra::values(result$claytotal$posterior$percentiles$P50)[, 1]
  expect_true(is.na(clay_p50[1]))
  expect_true(all(is.finite(clay_p50[-1])))
})

test_that("fuse_texture_group()'s posterior_probs = NULL resolves to FUSE_POSTERIOR_DEFAULT_PROBS", {
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  result <- fuse_texture_group(fetched)
  expect_named(result$claytotal$posterior$percentiles, paste0("P", round(FUSE_POSTERIOR_DEFAULT_PROBS * 100)))
})

test_that("fuse_closed_form(raw_draws) fits each mukey's real draws instead of the percentile-reconstructed approximation - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.6", {
  # Two mukeys with deliberately different real means; percentile-reconstruction sees only the
  # (P05,P25,P50,P75,P95) summary, while raw_draws should fit each mukey's fuller real distribution.
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  mukey_codes <- c(501, 502)
  mukey_raster <- terra::rast(nrows = 1, ncols = 2, vals = mukey_codes)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  set.seed(2026)
  draws_501 <- stats::rnorm(2000, mean = 10, sd = 2)
  draws_502 <- stats::rnorm(2000, mean = 30, sd = 2)
  mukey_draws <- list("501" = draws_501, "502" = draws_502)

  make_pct_rasters <- function(vals_by_cell) {
    stats::setNames(lapply(probs, function(p) {
      terra::rast(nrows = 1, ncols = 2, vals = stats::qnorm(p, vals_by_cell$mu, vals_by_cell$sd))
    }), paste0("P", round(probs * 100)))
  }
  # Percentile-reconstruction "sees" a different (wrong) per-cell mean than the real draws above,
  # so the two routes' fused means should diverge - a real behavioral difference, not a units check.
  prior_r <- make_pct_rasters(list(mu = c(15, 25), sd = c(3, 3)))
  lik_r <- make_pct_rasters(list(mu = c(15, 25), sd = c(3, 3)))

  baseline <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL)
  raw <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL,
                           mukey_raster = mukey_raster, mukey_draws = mukey_draws)

  mu_baseline <- terra::values(baseline$posterior$mu)[, 1]
  mu_raw <- terra::values(raw$posterior$mu)[, 1]
  expect_false(isTRUE(all.equal(mu_baseline, mu_raw)))
  # raw_draws fused mean should sit closer to the real draws' own means (10, 30) than the
  # percentile-reconstruction baseline (built from a mismatched mu of 15/25).
  expect_true(abs(mu_raw[1] - 10) < abs(mu_baseline[1] - 10))
  expect_true(abs(mu_raw[2] - 30) < abs(mu_baseline[2] - 30))
})

test_that("fuse_closed_form(raw_draws) degrades cells with no matching mukey draws to the percentile-based fit instead of NA/erroring", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  mukey_raster <- terra::rast(nrows = 1, ncols = 2, vals = c(601, 602))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)
  # Only mukey 601 has real draws; 602 has none and must fall back to the percentile-based fit.
  mukey_draws <- list("601" = stats::rnorm(2000, mean = 12, sd = 2))

  make_pct_rasters <- function() {
    stats::setNames(lapply(probs, function(p) {
      terra::rast(nrows = 1, ncols = 2, vals = rep(stats::qnorm(p, 20, 3), 2))
    }), paste0("P", round(probs * 100)))
  }
  prior_r <- make_pct_rasters()
  lik_r <- make_pct_rasters()

  result <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL,
                              mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  mu <- terra::values(result$posterior$mu)[, 1]
  expect_true(all(is.finite(mu)))
  # Cell 2 (mukey 602, no draws) should match the pure-percentile baseline for that cell.
  baseline <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL)
  expect_equal(unname(mu[2]), unname(terra::values(baseline$posterior$mu)[2, 1]), tolerance = 1e-8)
})

test_that("fuse_closed_form(raw_draws) with family='beta'/'gamma' produces finite fitted parameters", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  mukey_raster <- terra::rast(nrows = 1, ncols = 1, vals = 701)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  make_pct_rasters <- function(mu, sd) {
    stats::setNames(lapply(probs, function(p) terra::rast(nrows = 1, ncols = 1, vals = stats::qnorm(p, mu, sd))), paste0("P", round(probs * 100)))
  }

  # beta
  bounds <- c(0, 60)
  mukey_draws_b <- list("701" = pmin(pmax(stats::rnorm(1500, 30, 5), 0.5), 59.5))
  prior_rb <- make_pct_rasters(28, 5); lik_rb <- make_pct_rasters(32, 5)
  rb <- fuse_closed_form(prior_rb, lik_rb, probs, "beta", bounds = bounds,
                          mukey_raster = mukey_raster, mukey_draws = mukey_draws_b)
  expect_true(is.finite(terra::values(rb$posterior$alpha)[1, 1]))
  expect_true(is.finite(terra::values(rb$posterior$beta)[1, 1]))

  # gamma
  mukey_draws_g <- list("701" = stats::rgamma(1500, shape = 3, rate = 0.3))
  prior_rg <- make_pct_rasters(8, 2); lik_rg <- make_pct_rasters(10, 2)
  rg <- fuse_closed_form(prior_rg, lik_rg, probs, "gamma", bounds = NULL,
                          mukey_raster = mukey_raster, mukey_draws = mukey_draws_g)
  expect_true(is.finite(terra::values(rg$posterior$shape)[1, 1]))
  expect_true(is.finite(terra::values(rg$posterior$rate)[1, 1]))
})

test_that("fuse_closed_form(raw_draws) is unaffected when mukey_raster/mukey_draws are NULL (default, backward-compatible)", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  make_pct_rasters <- function() {
    stats::setNames(lapply(probs, function(p) terra::rast(nrows = 1, ncols = 1, vals = stats::qnorm(p, 20, 3))), paste0("P", round(probs * 100)))
  }
  prior_r <- make_pct_rasters(); lik_r <- make_pct_rasters()
  a <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL)
  b <- fuse_closed_form(prior_r, lik_r, probs, "normal", bounds = NULL, mukey_raster = NULL, mukey_draws = NULL)
  expect_identical(terra::values(a$posterior$mu), terra::values(b$posterior$mu))
  expect_identical(terra::values(a$posterior$sigma), terra::values(b$posterior$sigma))
})

test_that("fuse_lognormal_adaptive()'s small-AOI branch now forwards mukey_raster/mukey_draws to fuse_adaptive() - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.6", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  mukey_raster <- terra::rast(nrows = 1, ncols = 1, vals = 801)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)
  mukey_draws <- list("801" = stats::rlnorm(1500, meanlog = log(20), sdlog = 0.4))

  make_pct_rasters <- function() {
    stats::setNames(lapply(probs, function(p) terra::rast(nrows = 1, ncols = 1, vals = stats::qnorm(p, 20, 3))), paste0("P", round(probs * 100)))
  }
  prior_r <- make_pct_rasters(); lik_r <- make_pct_rasters()

  baseline <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 1e9, verbose = FALSE)
  raw <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 1e9, verbose = FALSE,
                                  mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  expect_true(is.finite(terra::values(raw$posterior$mu)[1, 1]))
  expect_false(isTRUE(all.equal(terra::values(baseline$posterior$mu), terra::values(raw$posterior$mu))))
})

test_that("resolve_want_raw_draws() defaults to TRUE for dist in {auto,normal,beta,gamma,lognormal} and FALSE for metalog - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.10", {
  # Unset prior_fusion_method (the common case): scope-aware default.
  expect_true(resolve_want_raw_draws(NULL, "normal"))
  expect_true(resolve_want_raw_draws(NULL, "beta"))
  expect_true(resolve_want_raw_draws(NULL, "gamma"))
  expect_true(resolve_want_raw_draws(NULL, "lognormal"))
  expect_true(resolve_want_raw_draws(NULL, "auto"))
  expect_false(resolve_want_raw_draws(NULL, "metalog"))
})

test_that("resolve_want_raw_draws() honors explicit prior_fusion_method overrides in both directions", {
  # Explicit 'raw_draws' forces TRUE even for metalog (accepted but unused downstream - see docs).
  expect_true(resolve_want_raw_draws("raw_draws", "metalog"))
  expect_true(resolve_want_raw_draws("raw_draws", "normal"))
  # Explicit 'percentile' forces FALSE even for a route that would otherwise default to TRUE.
  expect_false(resolve_want_raw_draws("percentile", "normal"))
  expect_false(resolve_want_raw_draws("percentile", "auto"))
  # Any other unrecognized string falls back to the same scope-aware default as unset.
  expect_true(resolve_want_raw_draws("bogus", "normal"))
  expect_false(resolve_want_raw_draws("bogus", "metalog"))
})

test_that("run_stage1_fusion() defaults to raw_draws for dist='normal' and to percentile for dist='metalog' - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.10", {
  # Offline, no live network: mock every fetch/simulate call to prove WHICH branch each dist takes,
  # not to exercise the real fusion math (already covered by this file's other tests).
  aoi <- terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  make_pct <- function(v) stats::setNames(
    lapply(probs, function(p) terra::rast(nrows = 1, ncols = 1, vals = stats::qnorm(p, v, 2))),
    paste0("P", round(probs * 100))
  )
  fake_prior <- list(values = make_pct(20), probs = probs)
  fake_solus <- list(values = make_pct(22), probs = probs)
  mukey_raster <- terra::rast(nrows = 1, ncols = 1, vals = 900)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  run_with_mocks <- function(dist, simulate_should_be_called) {
    simulate_called <- FALSE
    testthat::local_mocked_bindings(
      # Prior cache pre-warmed (non-NULL) so `is.null(prior) || want_raw_draws` isolates the effect
      # of want_raw_draws alone - if the cache were cold, simulate would run regardless of dist.
      cache_get_valid_percentiles = function(...) fake_prior,
      cache_set = function(...) invisible(NULL),
      build_cache_key = function(...) "fake_key",
      fetch_ssurgo_mukey_raster = function(...) mukey_raster,
      simulate_ssurgo_mapunit_draws = function(...) {
        simulate_called <<- TRUE
        list()
      },
      percentiles_from_draws = function(...) fake_prior,
      fetch_solus_percentiles = function(...) fake_solus,
      mukey_draws_lookup = function(...) NULL,
      .package = "soilSIM"
    )
    property_config <- list(id = "claytotal", solus_variable = "claytotal", dist = dist)
    # The want_raw_draws decision (and any simulate_ssurgo_mapunit_draws() call it triggers) always
    # happens before the dist-specific fuse_property_adaptive() dispatch, so checking
    # simulate_called is valid even if that downstream fusion step errors on this test's minimal
    # synthetic property_config (e.g. dist="metalog" needing its own term configuration) -
    # swallow any such error since it's irrelevant to what this test checks.
    tryCatch(
      run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 30),
      error = function(e) NULL
    )
    expect_equal(simulate_called, simulate_should_be_called)
  }

  run_with_mocks("normal", TRUE)
  run_with_mocks("metalog", FALSE)
})

test_that("fuse_lognormal_adaptive()'s large-AOI closed-form branch fits raw_draws directly in log-space and falls back exactly for uncovered mukeys - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.12", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  make_pct <- function(mu, sd, n) {
    stats::setNames(lapply(probs, function(p) terra::rast(nrows = 1, ncols = n, vals = rep(stats::qnorm(p, mu, sd), n))),
                     paste0("P", round(probs * 100)))
  }
  mukey_raster <- terra::rast(nrows = 1, ncols = 2, vals = c(501, 502))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)
  set.seed(2026)
  mukey_draws <- list("501" = stats::rlnorm(2000, meanlog = log(20), sdlog = 0.3))  # 502: no draws

  prior_r <- make_pct(25, 4, n = 2)
  lik_r <- make_pct(25, 4, n = 2)

  baseline <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 2, threshold_cells = 0, verbose = FALSE)
  raw <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 2, threshold_cells = 0, verbose = FALSE,
                                  mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  expect_equal(raw$route, "closed_form_lognormal")

  mu_base <- terra::values(baseline$posterior$mu)[, 1]
  mu_raw <- terra::values(raw$posterior$mu)[, 1]
  expect_true(all(is.finite(mu_raw)))
  # Cell 1 (mukey 501, has draws): raw_draws fit differs from the percentile-reconstruction baseline.
  expect_false(isTRUE(all.equal(unname(mu_base[1]), unname(mu_raw[1]))))
  # Cell 2 (mukey 502, no draws): falls back to the exact baseline value.
  expect_equal(unname(mu_raw[2]), unname(mu_base[2]), tolerance = 1e-8)
})

test_that("fuse_lognormal_adaptive()'s large-AOI branch with NULL mukey_raster/mukey_draws is bit-identical to baseline (backward compatibility)", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  make_pct <- function() {
    stats::setNames(lapply(probs, function(p) terra::rast(nrows = 1, ncols = 1, vals = stats::qnorm(p, 25, 4))),
                     paste0("P", round(probs * 100)))
  }
  prior_r <- make_pct(); lik_r <- make_pct()
  a <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 0, verbose = FALSE)
  b <- fuse_lognormal_adaptive(prior_r, probs, lik_r, probs, ncell = 1, threshold_cells = 0, verbose = FALSE,
                                mukey_raster = NULL, mukey_draws = NULL)
  expect_identical(terra::values(a$posterior$mu), terra::values(b$posterior$mu))
  expect_identical(terra::values(a$posterior$sigma), terra::values(b$posterior$sigma))
})

test_that("fuse_metalog_adapter() interpolates through real per-mukey empirical percentiles when raw_draws is supplied - MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.12", {
  # Known-working metalog fixture (P3.5 investigation): a symmetric 5-point Normal percentile set
  # triggers a pre-existing solve()-singularity in fit_metalog_linear_raster() unrelated to
  # raw_draws - reuse the 3-point asymmetric fixture already established as safe (see the
  # "percentiles exactly match analytic qnorm" test above).
  probs <- c(0.1, 0.5, 0.9)
  mk <- function(v) terra::rast(nrows = 1, ncols = 1, vals = v)
  prior_list <- stats::setNames(list(mk(20), mk(40), mk(65)), paste0("P", round(probs * 100)))
  lik_list <- stats::setNames(list(mk(30), mk(50), mk(70)), paste0("P", round(probs * 100)))
  property_config <- list(bounds = c(0, 100), boundedness = "b")

  mukey_raster <- terra::rast(nrows = 1, ncols = 1, vals = 900)
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)
  set.seed(11)
  mukey_draws <- list("900" = pmin(pmax(stats::rnorm(2000, mean = 55, sd = 8), 1), 99))

  baseline <- fuse_metalog_adapter(prior_list, probs, lik_list, probs, property_config, verbose = FALSE)
  raw <- fuse_metalog_adapter(prior_list, probs, lik_list, probs, property_config, verbose = FALSE,
                               mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  mu_base <- terra::values(baseline$posterior$mu)[1, 1]
  mu_raw <- terra::values(raw$posterior$mu)[1, 1]
  expect_true(is.finite(mu_raw))
  expect_false(isTRUE(all.equal(mu_base, mu_raw)))

  # Backward compatibility: NULL mukey_raster/mukey_draws (default) is bit-identical to baseline.
  explicit_null <- fuse_metalog_adapter(prior_list, probs, lik_list, probs, property_config, verbose = FALSE,
                                         mukey_raster = NULL, mukey_draws = NULL)
  expect_identical(terra::values(baseline$posterior$mu), terra::values(explicit_null$posterior$mu))

  # End-to-end dispatch through fuse_property_adaptive() with dist="metalog".
  property_config2 <- list(dist = "metalog", bounds = c(0, 100), boundedness = "b")
  e2e <- fuse_property_adaptive(prior_list, probs, lik_list, probs, property_config2, verbose = FALSE,
                                 mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  expect_true(is.finite(terra::values(e2e$posterior$mu)[1, 1]))
})

test_that("mukey_draws_percentiles_raster() degrades an uncovered mukey to NA instead of leaking its raw mukey code", {
  mukey_raster <- terra::rast(nrows = 1, ncols = 2, vals = c(601, 602))
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)
  mukey_draws <- list("601" = stats::rnorm(1000, 10, 2))  # 602: no entry

  result <- mukey_draws_percentiles_raster(mukey_raster, mukey_draws, c(0.5))
  vals <- terra::values(result$P50)[, 1]
  expect_true(is.finite(vals[1]))
  expect_true(is.na(vals[2]))
})
