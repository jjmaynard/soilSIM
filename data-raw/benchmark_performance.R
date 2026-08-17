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

if (identical(environment(), globalenv())) {
  benchmark_ssurgo_simulation()
  benchmark_texture_group_fusion()
  benchmark_monte_carlo_integration()
  benchmark_gp_hyperparameter_optimization()
  benchmark_fuse_general_kde()
  benchmark_simulate_cokey_generalized()
  benchmark_merge_adjusted_data()
  benchmark_apply_cross_property_constraints()
  benchmark_check_property_data_availability()
  cat("All benchmarks done.\n")
}
