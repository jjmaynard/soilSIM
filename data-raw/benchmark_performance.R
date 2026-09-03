# Maintainer-run performance benchmark harness - NOT part of the package build (data-raw/ is
# excluded via .Rbuildignore) and NOT part of the automated testthat suite (performance numbers
# are hardware-dependent and would make CI flaky). Run manually before/after a perf-related change
# to get a real before/after number, per PERFORMANCE_IMPROVEMENT_PLAN.md's methodology.
#
# Run from within a checkout of the soilSIM package directory (working directory = soilSIM/),
# e.g. via `Rscript data-raw/benchmark_performance.R` from inside soilSIM/. In this dev
# environment, always `unset PROJ_LIB` in the same shell invocation before running (a documented,
# terra-related environment quirk - see HANDOFF_NOTES.md/.onLoad()).
#
# Each benchmark function is self-contained (prints its own label + timing) so you can source this
# file and call just the one you need, e.g.:
#   source("data-raw/benchmark_performance.R"); benchmark_ssurgo_simulation()

devtools::load_all(".", quiet = TRUE)

# A ~200m x 200m box inside the same Salinas Valley area used elsewhere in this package's tests/
# vignettes, chosen for fast iteration (few mukeys/cokeys) while still exercising real live data.
small_wkt <- "POLYGON((-121.652 36.610, -121.650 36.610, -121.650 36.612, -121.652 36.612, -121.652 36.610))"

small_aoi <- function() {
  aoi <- terra::vect(small_wkt, crs = "epsg:4326")
  terra::project(aoi, "epsg:5070")
}

# ---------------------------------------------------------------------------
# 1. SSURGO tabular simulation pipeline (simulate_cokey_generalized(), GP depth-trend adjustment)
# ---------------------------------------------------------------------------
benchmark_ssurgo_simulation <- function() {
  cat("== [1] simulate_ssurgo_mapunit_draws() (small AOI) ==\n")
  aoi <- small_aoi()
  t <- system.time(draws <- simulate_ssurgo_mapunit_draws(aoi, top_depth = 0, bottom_depth = 5))
  print(t)
  cat("nrow(draws):", if (!is.null(draws)) nrow(draws) else "NULL", "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 2. Texture-group raster fusion (exercises fuse_texture_group()'s per-cell Monte Carlo path)
# ---------------------------------------------------------------------------
benchmark_texture_group_fusion <- function() {
  cat("== [2] run_stage1_fusion_group() for texture (small AOI) ==\n")
  aoi <- small_aoi()
  cache_dir <- tools::R_user_dir("soilSIM", "cache")
  if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE)

  composition_groups <- list(texture = list(members = c("clay", "sand", "silt")))
  property_configs <- list(
    clay = list(id = "clay", solus_variable = "claytotal", composition_group = "texture"),
    sand = list(id = "sand", solus_variable = "sandtotal", composition_group = "texture"),
    silt = list(id = "silt", solus_variable = "silttotal", composition_group = "texture")
  )
  t <- system.time(
    result <- run_stage1_fusion_group(aoi, "texture", composition_groups, property_configs,
                                       top_depth = 0, bottom_depth = 5)
  )
  print(t)
  cat("names(result):", if (!is.null(result)) paste(names(result), collapse = ", ") else "NULL", "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 3. Tabular Monte Carlo simulation (exercises apply_sum_constraints()) - uses the same cached
#    real Amador-area SSURGO data the vignettes use, so this runs offline/deterministically.
#
#    NOTE: `integrate_monte_carlo_with_gp()` (which reaches the Tier 1 `process_single_cokey()`
#    hot path via `process_cokeys_sequential()`/`process_cokeys_parallel()`) is deliberately not
#    exercised here - it has zero internal callers (confirmed via
#    `grep -rn "integrate_monte_carlo_with_gp(" R/`, only this doc cross-reference:
#    R/soilSIM-package.R:56) and its real usage contract (per
#    docs/05_gp_modeling_multivariate_adjustment.md) requires a `cokey_mapping` from
#    `match_soils_to_gp_models()` plus Monte Carlo realizations generated under the same
#    canonical property names (`clay_pct`/`sand_pct`/...) the GP models are trained on, not the
#    raw SSURGO names (`claytotal`/`sandtotal`/...) `generate_monte_carlo_realizations()` uses
#    here. When the `process_single_cokey()` fix (Tier 1) is implemented, benchmark/test it via a
#    purpose-built call using that real contract instead of stretching this generic harness to fit.
# ---------------------------------------------------------------------------
benchmark_monte_carlo_integration <- function() {
  cat("== [3] generate_monte_carlo_realizations() (cached Amador data) ==\n")
  cached_path <- system.file("extdata", "ssurgo_amador.rds", package = "soilSIM")
  if (!nzchar(cached_path)) {
    cat("Cached Amador data not found (package not installed?) - skipping.\n\n")
    return(invisible(NULL))
  }
  ssurgo_amador <- readRDS(cached_path)
  raw_data <- ssurgo_amador$ssurgo_data
  processed <- process_ssurgo_data(raw_data, max_depth = 150, verbose = FALSE)
  horizon_data <- processed$processed_data
  properties <- c("sandtotal", "claytotal", "silttotal", "dbovendry")
  infilled <- process_soil_properties_comprehensive(horizon_data, properties = properties, verbose = FALSE)

  t_mc <- system.time(
    mc_result <- generate_monte_carlo_realizations(
      soil_data = infilled, properties = properties, n_realizations = 200,
      simulation_config = list(max_depth = 150, auto_correlation = TRUE, correlation_fallback = "kssl_global"),
      parallel = FALSE, seed = 123
    )
  )
  cat("generate_monte_carlo_realizations():\n"); print(t_mc)
  cat("\n")
  invisible(list(mc = t_mc))
}

# ---------------------------------------------------------------------------
# 4. GP model training with hyperparameter optimization (exercises fit_individual_gp_model()/
#    optimize_gp_hyperparameters()'s ~16 GP_fit() calls per (property, group) pair)
# ---------------------------------------------------------------------------
benchmark_gp_hyperparameter_optimization <- function() {
  cat("== [4] build_stratified_gp_models(optimize_hyperparameters = TRUE) (cached Amador data) ==\n")
  cached_path <- system.file("extdata", "ssurgo_amador.rds", package = "soilSIM")
  if (!nzchar(cached_path)) {
    cat("Cached Amador data not found (package not installed?) - skipping.\n\n")
    return(invisible(NULL))
  }
  ssurgo_amador <- readRDS(cached_path)
  raw_data <- ssurgo_amador$ssurgo_data
  processed <- process_ssurgo_data(raw_data, max_depth = 150, verbose = FALSE)
  horizon_data <- processed$processed_data
  gp_train <- prepare_nrcs_training_data(horizon_data, max_depth = 150)

  t <- system.time(
    gp_models <- build_stratified_gp_models(
      gp_train, properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
      min_profiles_per_group = 3, min_observations_per_group = 15, optimize_hyperparameters = TRUE
    )
  )
  print(t)
  cat("total_models:", gp_models$model_summary$total_models, "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 5. General-KDE raster fusion (exercises fuse_general_kde()'s per-cell approxfun()+density()
#    path, reached via fuse_adaptive() whenever ncell <= threshold_cells). Synthetic percentile
#    rasters, bypassing network fetch - same style as the synthetic-raster fuse_texture_group()
#    benchmark referenced in PERFORMANCE_IMPROVEMENT_PLAN.md's Tier 1 log, not the
#    network-hitting benchmark_texture_group_fusion() above.
#
#    Added for PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4: the plan's own "Confirmed NOT bugs"
#    section previously waved this function off as "already gated behind threshold_cells -
#    working as designed" with no benchmark behind that call - the one unverified assumption in
#    an otherwise fully-benchmarked audit. This closes that gap.
# ---------------------------------------------------------------------------
make_synthetic_percentile_rasters <- function(ncell, base_p5, base_p50, base_p95, seed = 42) {
  nr <- max(1, floor(sqrt(ncell)))
  nc <- ceiling(ncell / nr)
  set.seed(seed)
  n <- nr * nc
  jitter <- function(v) v + stats::runif(n, -1, 1)
  list(
    P5  = terra::rast(nrows = nr, ncols = nc, vals = jitter(base_p5)),
    P50 = terra::rast(nrows = nr, ncols = nc, vals = jitter(base_p50)),
    P95 = terra::rast(nrows = nr, ncols = nc, vals = jitter(base_p95))
  )
}

benchmark_fuse_general_kde <- function(ncell = 2000) {
  prior_r <- make_synthetic_percentile_rasters(ncell, 20, 30, 45, seed = 42)
  lik_r   <- make_synthetic_percentile_rasters(ncell, 22, 32, 48, seed = 99)
  # nr*nc (floor(sqrt(ncell)) x ceiling(ncell/nr)) doesn't land exactly on the requested ncell -
  # force the general-KDE route off the ACTUAL raster cell count, not the requested one, or a
  # threshold_cells set from the requested value can silently fall short and route to
  # closed-form instead (confirmed happening with the naive `threshold_cells = ncell + 1` here).
  actual_ncell <- terra::ncell(prior_r[[1]])
  cat("== [5] fuse_general_kde() via fuse_adaptive() (", actual_ncell, "actual cells, synthetic) ==\n")
  t <- system.time(
    result <- fuse_adaptive(
      prior_r, lik_r, percentile_probs = c(0.05, 0.5, 0.95),
      family = "normal", threshold_cells = actual_ncell + 1, verbose = FALSE
    )
  )
  print(t)
  cat("route:", result$route, "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 6. simulate_cokey_generalized() (per-row Cholesky decomposition + correlated multivariate
#    draw + ILR texture transform) - a Tier 4 candidate, unlike most already-fixed row-loops its
#    per-row body is genuinely CPU-heavy, not just cheap dispatch overhead.
# ---------------------------------------------------------------------------
make_synthetic_property_correlation_matrices <- function(genhz = c("A", "B")) {
  vars <- c("db", "ph", "ilr1", "ilr2")
  m <- diag(4)
  dimnames(m) <- list(vars, vars)
  m["db", "ph"] <- m["ph", "db"] <- 0.3
  m["ilr1", "ilr2"] <- m["ilr2", "ilr1"] <- 0.2
  stats::setNames(lapply(genhz, function(g) m), genhz)
}

make_synthetic_texture_correlation_matrices <- function(genhz = c("A", "B")) {
  m <- matrix(c(1, -0.4, -0.4, -0.4, 1, -0.3, -0.4, -0.3, 1), nrow = 3)
  stats::setNames(lapply(genhz, function(g) m), genhz)
}

make_synthetic_sim_cokey_rows <- function(n_rows, seed = 1) {
  set.seed(seed)
  data.frame(
    genhz = sample(c("A", "B"), n_rows, replace = TRUE),
    compname = "benchseries", mukey = "1",
    cokey = as.character(rep(seq_len(ceiling(n_rows / 5)), each = 5)[seq_len(n_rows)]),
    hzdept_r = 0, hzdepb_r = 20,
    sim_comppct = sample(15:60, n_rows, replace = TRUE),
    dbovendry_l = 1.2, dbovendry_r = 1.4, dbovendry_h = 1.6,
    ph1to1h2o_l = 5.5, ph1to1h2o_r = 6.0, ph1to1h2o_h = 6.5,
    sandtotal_l = 30, sandtotal_r = 40, sandtotal_h = 50,
    silttotal_l = 25, silttotal_r = 35, silttotal_h = 45,
    claytotal_l = 15, claytotal_r = 25, claytotal_h = 35,
    stringsAsFactors = FALSE
  )
}

benchmark_simulate_cokey_generalized <- function(n_rows = 2000) {
  cat("== [6] simulate_cokey_generalized() (", n_rows, "synthetic horizon rows) ==\n")
  sim_cokey <- make_synthetic_sim_cokey_rows(n_rows)
  corr <- make_synthetic_property_correlation_matrices()
  txt_corr <- make_synthetic_texture_correlation_matrices()
  t <- system.time(result <- simulate_cokey_generalized(sim_cokey, corr, txt_corr))
  print(t)
  cat("nrow(result):", if (!is.null(result)) nrow(result) else "NULL", "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 7. merge_adjusted_data() (per-row which()-based full-table scan "join") - reached per-cokey via
#    apply_gp_depth_trends()/apply_nrcs_trend_adjustments(), same anti-pattern shape as the
#    already-fixed process_single_cokey() rescan (Tier 1, 11.4x).
# ---------------------------------------------------------------------------
make_synthetic_cokey_data <- function(n_depths, n_sims, seed = 4) {
  set.seed(seed)
  df <- expand.grid(hzdept_r = seq_len(n_depths), simulation_number = seq_len(n_sims))
  df$clay_pct <- stats::runif(nrow(df), 10, 40)
  df$sand_pct <- stats::runif(nrow(df), 20, 60)
  df
}

benchmark_merge_adjusted_data <- function(n_depths = 10, n_sims = 1000) {
  cat("== [7] merge_adjusted_data() (", n_depths, "depths x", n_sims, "realizations) ==\n")
  cokey_data <- make_synthetic_cokey_data(n_depths, n_sims)
  adjusted_data <- cokey_data
  adjusted_data$clay_pct <- adjusted_data$clay_pct + 1
  t <- system.time(
    result <- merge_adjusted_data(cokey_data, adjusted_data, c("clay_pct", "sand_pct"))
  )
  print(t)
  cat("nrow(result):", nrow(result), "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 8. apply_cross_property_constraints() (per-row texture sum/rescale) - same shape as the
#    already-fixed related_property_estimation() texture branch (Tier 2, 195x).
#
#    NOTE: correct_distribution_shapes() (this function's only caller) has zero internal call
#    sites anywhere in R/ (confirmed via `grep -rn "correct_distribution_shapes(" R/` - only a
#    vignette article doc reference) - exported-API-only, same status as the Tier 3
#    hz_quant_prob_mukey() item. Benchmarked/fixed for correctness/consistency, but expect no
#    measurable win on a real AOI pipeline run.
# ---------------------------------------------------------------------------
benchmark_apply_cross_property_constraints <- function(n_rows = 50000) {
  cat("== [8] apply_cross_property_constraints() (", n_rows, "synthetic rows; zero internal callers, see note) ==\n")
  set.seed(5)
  data <- data.frame(
    sandtotal = stats::runif(n_rows, 20, 60),
    claytotal = stats::runif(n_rows, 10, 40),
    silttotal = stats::runif(n_rows, 10, 40)
  )
  t <- system.time(
    result <- apply_cross_property_constraints(data, c("sandtotal", "claytotal", "silttotal"))
  )
  print(t)
  cat("nrow(result):", nrow(result), "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 9. check_property_data_availability() (nested per-row/per-property NA check, re-indexing the
#    same column on every (row, property) pair instead of once) - real caller in
#    generate_monte_carlo_realizations()'s pipeline (monte-carlo.R:1840).
# ---------------------------------------------------------------------------
benchmark_check_property_data_availability <- function(n_rows = 20000) {
  cat("== [9] check_property_data_availability() (", n_rows, "synthetic rows) ==\n")
  set.seed(6)
  properties <- c("claytotal", "sandtotal", "silttotal", "dbovendry", "ph1to1h2o")
  soil_data <- as.data.frame(stats::setNames(
    lapply(properties, function(p) {
      v <- stats::runif(n_rows, 0, 50)
      v[sample(n_rows, floor(n_rows * 0.1))] <- NA
      v
    }),
    paste0(properties, "_r")
  ))
  t <- system.time(
    result <- check_property_data_availability(soil_data, properties, config = list())
  )
  print(t)
  cat("sum(has_data):", sum(result), "\n\n")
  invisible(t)
}

# ---------------------------------------------------------------------------
# 10. fuse_general_kde()'s raw_draws route vs the default percentile-reconstruction route
#     (MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.9) - raw_draws does a per-cell mukey lookup
#     against real Monte Carlo draws instead of reconstructing from a 5-percentile summary, so
#     it's expected to cost more; this quantifies by how much, at a representative cell count and
#     mukey cardinality, using the same synthetic-raster style as benchmark #5 (no live network).
# ---------------------------------------------------------------------------
benchmark_fuse_general_kde_raw_draws <- function(ncell = 2000, n_unique_mukeys = 8, draws_per_mukey = 1500) {
  prior_r <- make_synthetic_percentile_rasters(ncell, 20, 30, 45, seed = 42)
  lik_r   <- make_synthetic_percentile_rasters(ncell, 22, 32, 48, seed = 99)
  actual_ncell <- terra::ncell(prior_r[[1]])

  # Assign each cell one of n_unique_mukeys codes (a realistic AOI has far fewer unique mukeys
  # than cells - mukey_prior_samples is precomputed once per unique mukey, not once per cell, so
  # mukey cardinality - not cell count - is what should drive any extra cost here).
  set.seed(7)
  mukey_codes <- sample(seq_len(n_unique_mukeys), actual_ncell, replace = TRUE)
  mukey_raster <- prior_r[[1]]
  terra::values(mukey_raster) <- mukey_codes
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  mukey_draws <- stats::setNames(
    lapply(seq_len(n_unique_mukeys), function(i) {
      set.seed(100 + i)
      stats::rlnorm(draws_per_mukey, meanlog = log(25 + i), sdlog = 0.3)
    }),
    as.character(seq_len(n_unique_mukeys))
  )

  cat("== [10] fuse_general_kde() default vs raw_draws (", actual_ncell, "cells,",
      n_unique_mukeys, "unique mukeys, synthetic) ==\n")

  t_default <- system.time(
    result_default <- fuse_adaptive(
      prior_r, lik_r, percentile_probs = c(0.05, 0.5, 0.95),
      family = "normal", threshold_cells = actual_ncell + 1, verbose = FALSE
    )
  )
  cat("default route:\n"); print(t_default)

  t_raw <- system.time(
    result_raw <- fuse_adaptive(
      prior_r, lik_r, percentile_probs = c(0.05, 0.5, 0.95),
      family = "normal", threshold_cells = actual_ncell + 1, verbose = FALSE,
      mukey_raster = mukey_raster, mukey_draws = mukey_draws
    )
  )
  cat("raw_draws route:\n"); print(t_raw)

  ratio <- unname(t_raw[["elapsed"]] / t_default[["elapsed"]])
  cat(sprintf("raw_draws / default elapsed ratio: %.2fx\n\n", ratio))
  invisible(list(default = t_default, raw_draws = t_raw, ratio = ratio))
}

# ---------------------------------------------------------------------------
# 11. fuse_texture_group() default vs raw_draws (RASTER_STATISTICS_INTEGRATION_PLAN.md-adjacent -
#    MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md task P2.11). Mirrors benchmark_fuse_general_kde_raw_
#    draws() above, but for the joint texture-group route (fuse_texture_group_batch_core()'s
#    per-cell Cholesky/MC sampler), which P2.9 never benchmarked - only the single-property
#    general-KDE route was measured there.
# ---------------------------------------------------------------------------
make_synthetic_texture_percentile_rasters <- function(ncell, seed = 42) {
  nr <- max(1, floor(sqrt(ncell)))
  nc <- ceiling(ncell / nr)
  n <- nr * nc
  set.seed(seed)
  clay <- pmax(pmin(stats::rnorm(n, 20, 3), 55), 1)
  sand <- pmax(pmin(45 - 0.6 * (clay - 20) + stats::rnorm(n, 0, 2), 90), 1)
  silt <- pmax(100 - clay - sand, 1)
  mk_layer <- function(v) terra::rast(nrows = nr, ncols = nc, vals = v)
  list(clay = clay, sand = sand, silt = silt, nr = nr, nc = nc,
       template = mk_layer(clay))
}

benchmark_fuse_texture_group_batch_raw_draws <- function(ncell = 2000, n_unique_mukeys = 8, draws_per_mukey = 1500) {
  base <- make_synthetic_texture_percentile_rasters(ncell, seed = 42)
  actual_ncell <- length(base$clay)
  probs <- c(0.05, 0.5, 0.95)

  jitter_rasters <- function(v, spread, seed) {
    set.seed(seed)
    list(
      lo = terra::rast(nrows = base$nr, ncols = base$nc, vals = v - spread + stats::runif(length(v), -0.5, 0.5)),
      p50 = terra::rast(nrows = base$nr, ncols = base$nc, vals = v + stats::runif(length(v), -0.5, 0.5)),
      hi = terra::rast(nrows = base$nr, ncols = base$nc, vals = v + spread + stats::runif(length(v), -0.5, 0.5))
    )
  }
  clay_prior <- jitter_rasters(base$clay, 5, 1); clay_lik <- jitter_rasters(base$clay + 2, 5, 2)
  sand_prior <- jitter_rasters(base$sand, 6, 3); sand_lik <- jitter_rasters(base$sand + 2, 6, 4)
  silt_prior <- jitter_rasters(base$silt, 4, 5); silt_lik <- jitter_rasters(base$silt - 4, 4, 6)

  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  # Same "far fewer unique mukeys than cells" realism convention as benchmark_fuse_general_kde_raw_
  # draws() above - mukey cardinality, not cell count, is what should drive raw_draws' extra cost.
  set.seed(7)
  mukey_codes <- sample(seq_len(n_unique_mukeys), actual_ncell, replace = TRUE)
  mukey_raster <- base$template
  terra::values(mukey_raster) <- mukey_codes
  names(mukey_raster) <- "mukey"
  mukey_raster <- terra::as.factor(mukey_raster)

  mukey_texture_draws <- stats::setNames(
    lapply(seq_len(n_unique_mukeys), function(i) {
      set.seed(200 + i)
      n <- draws_per_mukey
      c_d <- pmax(stats::rnorm(n, 20 + i, 3), 1)
      s_d <- pmax(45 - 0.6 * (c_d - 20) + stats::rnorm(n, 0, 2), 1)
      si_d <- pmax(100 - c_d - s_d, 1)
      cbind(clay_total = c_d, sand_total = s_d, silt_total = si_d)
    }),
    as.character(seq_len(n_unique_mukeys))
  )

  cat("== [11] fuse_texture_group() default vs raw_draws (", actual_ncell, "cells,",
      n_unique_mukeys, "unique mukeys, synthetic) ==\n")

  t_default <- system.time(result_default <- fuse_texture_group(fetched))
  cat("default route:\n"); print(t_default)

  t_raw <- system.time(
    result_raw <- fuse_texture_group(fetched, mukey_raster = mukey_raster, mukey_texture_draws = mukey_texture_draws)
  )
  cat("raw_draws route:\n"); print(t_raw)

  ratio <- unname(t_raw[["elapsed"]] / t_default[["elapsed"]])
  cat(sprintf("raw_draws / default elapsed ratio: %.2fx\n\n", ratio))
  invisible(list(default = t_default, raw_draws = t_raw, ratio = ratio))
}

# ---------------------------------------------------------------------------
# A.6 - per-pixel ensemble re-marginalization (raster-fusion-bridge.R). Live: SDA + SOLUS.
#   `benchmark_perpixel_remarginalization()` - times the ensemble build, per-window fusion of the
#     Saxton-Rawls inputs (sand/silt/clay/db - SOLUS100 has NO water-retention variable), and
#     remarginalized_awc() whole-grid vs tiled (asserts they agree).
#   `compare_perpixel_vs_zonal(prop_id, solus_var)` - the "does sub-mukey spatial structure
#     matter?" number: within-mukey SD of the per-pixel re-marginalized P50 vs. between-mukey SD.
#     ratio << 1 => a mukey-collapsed (zonal) feed loses little => promote
#     `zonal_distribution_from_posterior()`. Default property: claytotal.
# ---------------------------------------------------------------------------
.awc_windows <- list(c(0, 5), c(5, 15), c(15, 30))

# Per-window run_stage1_fusion() posteriors for a set of {sim-column id -> SOLUS variable} configs.
.perpixel_posteriors <- function(aoi, windows, cfgs) {
  out <- list()
  for (nm in names(cfgs)) {
    per_w <- list()
    for (w in windows) {
      r <- tryCatch(run_stage1_fusion(aoi, list(id = nm, solus_variable = cfgs[[nm]], dist = "auto"),
                                      top_depth = w[1], bottom_depth = w[2]),
                    error = function(e) NULL)
      if (!is.null(r)) per_w[[paste0(w[1], "-", w[2])]] <- list(percentiles = r$posterior$percentiles)
    }
    if (length(per_w) > 0) out[[nm]] <- per_w
  }
  out
}

# Saxton-Rawls AWC inputs: sim-column name -> fetchSOLUS variable.
.awc_cfgs <- c(sand_total = "sandtotal", silt_total = "silttotal",
               clay_total = "claytotal", db = "dbovendry")

benchmark_perpixel_remarginalization <- function() {
  cat("== [A.6] extract_mukey_joint_ensemble() + run_stage1_fusion() + remarginalized_awc() ==\n")
  aoi <- small_aoi()
  windows <- .awc_windows

  t_ens <- system.time(ens <- extract_mukey_joint_ensemble(aoi, windows))
  cat("extract_mukey_joint_ensemble():\n"); print(t_ens)
  if (is.null(ens)) { cat("ensemble NULL (live data unavailable) - skipping.\n\n"); return(invisible(NULL)) }
  cat("ensemble properties present:", paste(ens$properties, collapse = ", "), "\n")

  t_post <- system.time(pbpw <- .perpixel_posteriors(aoi, windows, .awc_cfgs))
  cat("run_stage1_fusion() x (", length(.awc_cfgs), "props x", length(windows), "windows):\n"); print(t_post)
  if (!all(names(.awc_cfgs) %in% names(pbpw))) {
    cat("SOLUS fusion unavailable for", paste(setdiff(names(.awc_cfgs), names(pbpw)), collapse = ", "),
        "- skipping.\n\n"); return(invisible(NULL))
  }

  t_whole <- system.time(res_whole <- remarginalized_awc(ens, pbpw, n_out = 250))
  cat("remarginalized_awc() whole-grid (n_kept =", res_whole$n_kept, "):\n"); print(t_whole)
  t_tiled <- system.time(res_tiled <- remarginalized_awc(ens, pbpw, n_out = 250, tile_rows = 64))
  cat("remarginalized_awc() tiled (n_tiles =", res_tiled$n_tiles, "):\n"); print(t_tiled)
  d <- as.numeric(terra::global(abs(res_whole$awc_cm$P50 - res_tiled$awc_cm$P50), "max", na.rm = TRUE))
  cat(sprintf("max |whole - tiled| P50 AWC: %.3g (should be ~0)\n", d))
  m <- as.numeric(terra::global(res_whole$awc_cm$P50, "mean", na.rm = TRUE))
  cat(sprintf("mean P50 AWC over 0-30 cm: %.2f cm\n\n", m))
  invisible(list(ens = t_ens, post = t_post, whole = t_whole, tiled = t_tiled))
}

compare_perpixel_vs_zonal <- function(prop_id = "clay_total", solus_var = "claytotal",
                                      windows = .awc_windows) {
  cat(sprintf("== [A.6] within-mukey vs between-mukey spread of per-pixel %s (P50) ==\n", prop_id))
  aoi <- small_aoi()
  ens <- extract_mukey_joint_ensemble(aoi, windows)
  if (is.null(ens)) { cat("ensemble NULL - skipping.\n\n"); return(invisible(NULL)) }
  pbpw <- .perpixel_posteriors(aoi, windows, stats::setNames(solus_var, prop_id))
  if (is.null(pbpw[[prop_id]])) { cat("no fused posterior for", prop_id, "- skipping.\n\n"); return(invisible(NULL)) }

  rm_res <- remarginalize_ensemble_to_posterior(ens, pbpw, n_out = 250, summarize = TRUE,
                                                probs = c(0.5))
  ratios <- c()
  for (w in names(rm_res$percentiles[[prop_id]])) {
    p50 <- rm_res$percentiles[[prop_id]][[w]]$P50
    mk <- terra::resample(ens$mukey_raster, p50, method = "near")
    within  <- terra::zonal(p50, mk, fun = "sd",   na.rm = TRUE)
    between <- terra::zonal(p50, mk, fun = "mean", na.rm = TRUE)
    mw <- mean(within[[2]], na.rm = TRUE); sb <- stats::sd(between[[2]], na.rm = TRUE)
    ratios[w] <- mw / sb
    cat(sprintf("  window %-7s  mean within-mukey SD %.3f | between-mukey SD %.3f | ratio %.2f\n",
                w, mw, sb, mw / sb))
  }
  cat(sprintf("=> mean ratio across windows: %.2f  (<< 1 => zonal/mukey-collapsed feed loses little)\n\n",
              mean(ratios, na.rm = TRUE)))
  invisible(ratios)
}

# ---------------------------------------------------------------------------
# B.5 - analyze_soil_statistics() `_safe` chain vs the (now fixed) enhanced chain, on the same
#   cached Amador SSURGO data as benchmark #3. `fitdistrplus::fitdist()` parametric fitting cost
#   in the enhanced chain is undocumented and this measures it.
# ---------------------------------------------------------------------------
benchmark_statistics_chains <- function() {
  cat("== [B.5] analyze_soil_statistics(): _safe chain vs enhanced chain ==\n")
  cached_path <- system.file("extdata", "ssurgo_amador.rds", package = "soilSIM")
  if (!nzchar(cached_path)) { cat("Cached Amador data not found - skipping.\n\n"); return(invisible(NULL)) }

  ssurgo_amador <- readRDS(cached_path)
  processed <- process_ssurgo_data(ssurgo_amador$ssurgo_data, max_depth = 150, verbose = FALSE)$processed_data

  t_safe <- system.time(r_safe <- analyze_soil_statistics(processed, use_enhanced_chain = FALSE, verbose = FALSE))
  cat("_safe chain:\n"); print(t_safe)

  t_enh <- system.time(r_enh <- analyze_soil_statistics(processed, use_enhanced_chain = TRUE, verbose = FALSE))
  cat("enhanced chain:\n"); print(t_enh)

  ratio <- unname(t_enh[["elapsed"]] / max(t_safe[["elapsed"]], 1e-6))
  cat(sprintf("enhanced / _safe elapsed ratio: %.2fx | properties: %d/%d fitted\n\n",
              ratio,
              length(r_enh$distribution_analysis$fitted_distributions %||% list()),
              length(r_safe$distribution_analysis$fitted_distributions %||% list())))
  invisible(list(safe = t_safe, enhanced = t_enh, ratio = ratio))
}

if (identical(environment(), globalenv())) {
  benchmark_ssurgo_simulation()
  benchmark_texture_group_fusion()
  benchmark_monte_carlo_integration()
  benchmark_gp_hyperparameter_optimization()
  benchmark_fuse_general_kde()
  benchmark_fuse_general_kde_raw_draws()
  benchmark_simulate_cokey_generalized()
  benchmark_merge_adjusted_data()
  benchmark_apply_cross_property_constraints()
  benchmark_check_property_data_availability()
  benchmark_fuse_texture_group_batch_raw_draws()
  benchmark_perpixel_remarginalization()
  compare_perpixel_vs_zonal()
  benchmark_statistics_chains()
  cat("All benchmarks done.\n")
}
