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

if (identical(environment(), globalenv())) {
  benchmark_ssurgo_simulation()
  benchmark_texture_group_fusion()
  benchmark_monte_carlo_integration()
  benchmark_gp_hyperparameter_optimization()
  cat("All benchmarks done.\n")
}
