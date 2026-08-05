# GP Depth Modeling & Multivariate Adjustment

## Overview

This functional area covers two tightly coupled files. `gp-modeling.R`
prepares NRCS/SSURGO component-horizon data and fits Gaussian Process
(GP) models of how each soil property trends with depth, stratified by
soil group (series, taxonomic class, particle-size class, etc., with
hierarchical fallback to coarser groupings when a stratum is too
sparse). `multivariate-adjustment.R` takes the output of Monte Carlo
simulation (see the [Monte Carlo
Simulation](https://jjmaynard.github.io/soilSIM/articles/architecture-monte-carlo-simulation.md)
engine module) across an entire dataset and integrates it with these
fitted GP depth trends: for every simulated cokey it nudges each
realization at each depth toward the GP-predicted trend using a
reference-quantile (“rank-preserving”) transform, so that the depth
trend becomes realistic while the *within-depth* correlation structure
between properties (e.g. clay vs. bulk density vs. pH) is preserved
rather than scrambled by adjusting each property independently. The
result is a GP-depth-adjusted, correlation-preserving set of property
realizations per component (cokey), ready for downstream statistical
summarization and validation.

## Core Functions

### gp-modeling.R

#### `prepare_nrcs_training_data()`

``` r

prepare_nrcs_training_data(
  nrcs_combined_data,
  grouping_strategy = "auto",
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  target_min_groups = 3,
  max_depth = 250,
  validation_config = NULL
)
```

**Purpose**: Turn raw combined NRCS component/horizon data into a
cleaned, grouped training set ready for
[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md).
**Parameters**: - `nrcs_combined_data` - combined NRCS data with horizon
and component information (must contain `cokey`, `hzdept_r`). -
`grouping_strategy` - `"auto"` triggers
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md);
otherwise one of `"soil_series"`, `"taxonomic_class"`,
`"particle_size"`, `"soil_grtgroup"`, `"soil_suborder"`, `"soil_order"`,
`"none"`. - `min_profiles_per_group` / `min_observations_per_group` -
minimum distinct `cokey`s / total rows a group must have to be
considered adequate. - `target_min_groups` - minimum number of adequate
groups the auto-selector wants to find before accepting a strategy. -
`max_depth` - horizons with `hzdepb_r` beyond this (cm) are dropped. -
`validation_config` - quality-threshold config; defaults to
`get_default_configuration("validation")`. **Returns**: A data frame
with standardized property columns (`clay_pct`, `sand_pct`, `silt_pct`,
`pH`, `organic_matter`, `bulk_density`, `cec`, `awc`), a `soil_group`
column, `depth_midpoint`, and only rows belonging to groups meeting
`min_observations_per_group`. **Algorithm**: Runs
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
standardizes property names via
[`standardize_property_names()`](https://jjmaynard.github.io/soilSIM/reference/standardize_property_names.md),
computes `unsuitable_horizon` with
`is_unsuitable(processed_data, strict_mode = TRUE)` and filters those
horizons out (a `|>`-pipe-safe rewrite of the original
`%>% mutate(is_unsuitable(.))` pattern - see the inline code comment
about why the flag must be precomputed rather than referenced via a bare
`.` inside `mutate()`), auto-selects a grouping strategy if requested,
reconciles `hzdept_r`/`hzdepb_r` from whichever depth columns are
populated, coalesces synonym columns into the eight standard property
names via
[`safe_coalesce()`](https://jjmaynard.github.io/soilSIM/reference/safe_coalesce.md),
assigns `soil_group` via `create_soil_groups()`, drops inadequate
groups, and logs a processing summary.

#### `select_optimal_grouping()`

``` r

select_optimal_grouping(data, min_profiles, min_obs, target_groups)
```

**Purpose**: Pick the finest available soil-grouping strategy (series
before taxonomic class before particle size, etc.) that yields at least
`target_groups` adequate groups. **Parameters**: `data` - input NRCS
data; `min_profiles`/`min_obs` - per-group thresholds; `target_groups` -
minimum adequate-group count sought. **Returns**: The selected strategy
name (character), or the best-available strategy if none meets the
target. **Algorithm**: Iterates the hierarchy
`soil_series -> taxonomic_class -> particle_size -> soil_grtgroup -> soil_suborder -> soil_order -> none`,
skipping strategies whose required column is absent, builds a test
grouping via `create_test_grouping()`, computes per-group
profile/observation/ depth-range statistics, and returns the first
strategy meeting `n_adequate_groups >= target_groups`; falls back to
whichever strategy produced the most adequate groups.

#### `validate_training_groups()`

``` r

validate_training_groups(
  processed_data,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  min_profiles = 3,
  min_observations = 15
)
```

**Purpose**: Assess whether the grouped training data is adequate for
reliable GP model fitting. **Parameters**: `processed_data` - output of
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md);
`properties` - which columns to check coverage for;
`min_profiles`/`min_observations` - per-group adequacy thresholds.
**Returns**: A list with `overall_adequacy` (logical),
`group_validation` (per-group adequacy data frame), `property_coverage`
(per-property, per-group coverage data frames), and `recommendations`
(character vector of warnings, e.g. “Insufficient adequate groups” or
low-coverage notices). **Algorithm**: Summarizes `n_profiles`,
`n_observations`, and `depth_range` per `soil_group`, flags groups
meeting all three adequacy criteria, and for each requested property
computes per-group non-NA coverage rate among adequate groups, emitting
a recommendation when overall coverage falls below 30%.

#### `build_stratified_gp_models()` - Master Function

``` r

build_stratified_gp_models(
  processed_nrcs_data,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  optimize_hyperparameters = TRUE
)
```

**Purpose**: Build one GP depth-trend model per (property, soil-group)
combination, with automatic hierarchical fallback for undersized groups.
**Parameters**: `processed_nrcs_data` - output of
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md);
`properties` - properties to model;
`min_profiles_per_group`/`min_observations_per_group` - adequacy
thresholds used by the fallback grouping; `optimize_hyperparameters` -
whether to run cross-validated correlation-family selection
([`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md))
or fit a single default GP. **Returns**: A named list keyed by property,
each element
`list(type = "stratified_grouped", models = <list of per-group model objects>, grouping_summary, model_diagnostics, property, n_groups, total_training_points)`,
plus a `model_summary` element from `create_model_summary()`.
**Algorithm**: For each property, calls `apply_hierarchical_grouping()`
to assign every row a `final_group` (falling back through particle size
-\> great group -\> suborder -\> order -\>
`fallback_general`/`general_pool` when the original `soil_group` is too
small), summarizes each final group, and calls
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)
per group inside a `tryCatch` (failures are logged via
[`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md)
and skipped rather than aborting the whole property).

#### `fit_individual_gp_model()`

``` r

fit_individual_gp_model(data, property, optimize_hyperparameters = TRUE)
```

**Purpose**: Fit a single 1-D GP regression of a property against depth
for one soil group. **Parameters**: `data` - rows for one group;
`property` - column name to model; `optimize_hyperparameters` - use
[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)
(CV-selected correlation family) if `TRUE`, else a plain
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html).
**Returns**: `NULL` if fitting is not possible (fewer than 3 aggregated
depth points, no variance in values, or no depth range), otherwise
`list(model = list(gp_model, depth_scaling = list(min, max, range), training_data, n_training_points, property), diagnostics)`.
**Algorithm**: Aggregates to one mean value per `cokey`/depth then per
depth (`mean_value`, `n_profiles`, `sd_value`), filters degenerate rows,
min-max scales depths to `[0,1]` (`depth_scaling` is stored so
predictions can rescale new depths later), fits the GP via
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)
(optionally via the CV wrapper), and computes diagnostics with
`calculate_model_diagnostics()`.

#### `adjust_multivariate_depthwise_GP()`

``` r

adjust_multivariate_depthwise_GP(simulated_list, gp_models, depths, primary_property = NULL)
```

**Purpose**: The original (pre-refactor) multivariate
correlation-preserving adjustment routine, operating directly on
matrices rather than long-format simulation data frames. **Parameters**:
`simulated_list` - named list of `[depths x simulations]` matrices, one
per property; `gp_models` - named list of fitted GP models (or nested
`list(gp_model = ...)` structures) per property; `depths` - depth vector
matching matrix rows; `primary_property` - reference property whose
surface-depth quantile ordering is applied to every other property
(defaults to the first property in `simulated_list`). **Returns**: A
named list of adjusted matrices, same shape as `simulated_list`.
**Algorithm**: Validates that all matrices share dimensions and that
every property has a corresponding GP model, predicts GP means at each
depth (via
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md)
for nested models or
[`GPfit::predict.GP()`](https://rdrr.io/pkg/GPfit/man/predict.html)
directly, falling back to
[`rowMeans()`](https://rdrr.io/r/base/colSums.html) on prediction
failure), computes an ECDF of the primary property’s surface-depth
values to fix a per-simulation reference quantile, then depth-by-depth
(from depth 2 onward) computes a GP ratio between consecutive depth
means, nudges each simulation toward `previous_value * gp_ratio` at the
*same reference quantile* in the current property’s own distribution,
and re-maps the nudged values back onto the current property’s original
ECDF to keep its marginal distribution shape intact. Logs a post-hoc
check of maximum correlation-matrix difference (original vs. adjusted)
at the first five depths, warning if it exceeds 0.1. This is the
array-based sibling of `multivariate-adjustment.R`’s
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md),
which operates on the same algorithm but is invoked from the long-format
integration pipeline.

#### `validate_correlation_preservation()`

``` r

validate_correlation_preservation(original_list, adjusted_list, depths)
```

**Purpose**: Quantify how much the depth-wise correlation matrix between
properties changed between `original_list` and `adjusted_list` (the
matrix-format counterparts used by
[`adjust_multivariate_depthwise_GP()`](https://jjmaynard.github.io/soilSIM/reference/adjust_multivariate_depthwise_GP.md)).
**Returns**: A data frame with one row per depth (`depth_idx`,
`depth_value`, `max_cor_difference`, `mean_cor_difference`, the latter
averaged over the upper triangle of the correlation-difference matrix).
**Algorithm**: For each depth, extracts a `[simulations x properties]`
slice from both lists, computes `cor(..., use = "complete.obs")` for
each, and records the max/mean absolute difference; wrapped in per-depth
`tryCatch` so one bad depth does not abort the whole validation.

#### `simulate_soil_properties()`

``` r

simulate_soil_properties(
  target_cokey,
  nrcs_gp_models = NULL,
  cokey_mapping = NULL,
  sim_data,
  correlation_matrices,
  txt_correlation_matrices,
  n_simulations = 100,
  use_nrcs_gp = !is.null(nrcs_gp_models),
  preserve_correlations = TRUE
)
```

**Purpose**: End-to-end single-cokey convenience wrapper: generate Monte
Carlo realizations for one component, then apply either NRCS GP or local
GP depth-trend adjustment with correlation preservation, in one call.
**Parameters**: `target_cokey` - the cokey to simulate;
`nrcs_gp_models`/`cokey_mapping` - fitted NRCS GP models and their
cokey-to-group mapping (optional; enables NRCS-GP mode); `sim_data` -
processed per-horizon simulation input for this cokey (and others -
filtered internally);
`correlation_matrices`/`txt_correlation_matrices` - genhz-keyed (or
flat) correlation matrices; `n_simulations` - number of Monte Carlo
realizations; `use_nrcs_gp` - use NRCS GP models (defaults to `TRUE` iff
`nrcs_gp_models` supplied) vs. locally-fit GP models;
`preserve_correlations` - whether to run the correlation-preserving
adjustment path. **Returns**: A long-format data frame of enhanced
simulations for `target_cokey`, with `validation`, `n_simulations`,
`use_nrcs_gp`, and `preserve_correlations` attached as attributes (empty
data frame if the cokey has no data). **Algorithm**: Filters `sim_data`
to `target_cokey`; picks a single flat correlation matrix via
[`select_simulation_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/select_simulation_correlation_matrix.md)
(using the cokey’s modal `genhz` group, falling back to
`txt_correlation_matrices` or to empirically-estimated correlations);
calls
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)
and flattens its `[horizon, property, realization]` array to long format
via
[`flatten_simulation_array_to_long()`](https://jjmaynard.github.io/soilSIM/reference/flatten_simulation_array_to_long.md);
then, if NRCS GP models and a cokey mapping are supplied, looks up the
cokey’s `gp_model_group` (falling back to `"fallback_general"`) and
calls
[`apply_nrcs_gp_adjustments_with_correlations()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_gp_adjustments_with_correlations.md),
otherwise calls
[`apply_local_gp_adjustments_with_correlations()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments_with_correlations.md);
finally runs `validate_enhanced_simulations()` and attaches the result
as an attribute.

#### `match_soils_to_gp_models()`

``` r

match_soils_to_gp_models(
  simulated_cokeys,
  nrcs_combined_data,
  gp_models,
  property = "clay_pct",
  matching_strategy = "exact_cokey"
)
```

**Purpose**: Map every simulated cokey to the GP model group (from
[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md))
that should be used to adjust its depth trend. **Parameters**:
`simulated_cokeys` - cokeys needing a mapping; `nrcs_combined_data` -
source NRCS data carrying taxonomic columns; `gp_models` - fitted models
keyed by property; `property` - which property’s model groups to match
against; `matching_strategy` - currently only `"exact_cokey"` is
implemented (hierarchical taxonomic matching). **Returns**: A two-column
data frame, `sim_cokey` and `gp_model_group` (`"pooled_model"`
short-circuit if `gp_models[[property]]$type` is not
`"stratified_grouped"`). **Algorithm**: Builds a distinct lookup of
taxonomic columns (`compname`, `taxclname`, `taxpartsize`,
`taxgrtgroup`, `taxsuborder`, `taxorder`) from `nrcs_combined_data`,
left-joins it onto the target cokeys, then calls
`apply_matching_hierarchy()` to fill `gp_model_group` in priority order
(`compname` exact match -\> `taxclname` exact match -\>
fallback-prefixed particle size/great group/suborder/order matches -\>
`"fallback_general"` -\> `"general_pool"` -\> first available group as
last resort). Logs a per-group match-count summary and warns about any
still-unmatched cokeys.

#### `predict_gp_depth_trends()`

``` r

predict_gp_depth_trends(gp_model_info, new_depths)
```

**Purpose**: Predict a fitted GP model’s property values at arbitrary
new depths. **Parameters**: `gp_model_info` - the `model` element
produced by
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)
(contains `gp_model` and `depth_scaling`); `new_depths` - numeric vector
of depths to predict at. **Returns**: Numeric vector of predictions
(`Y_hat`), same length as `new_depths`; `NA`-filled if the model info is
invalid or prediction errors. **Algorithm**: Rescales `new_depths` into
the `[0,1]` range the model was trained on using the stored
`depth_scaling$min`/`range`, clamps to `[0,1]`, and calls
[`GPfit::predict.GP()`](https://rdrr.io/pkg/GPfit/man/predict.html).

#### `validate_gp_models()`

``` r

validate_gp_models(gp_models, nrcs_data, validation_depths = seq(0, 200, by = 10))
```

**Purpose**: Diagnostic pass over an entire fitted GP model set,
checking predictability, trend monotonicity, and realistic value ranges.
**Parameters**: `gp_models` - output of
[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md);
`nrcs_data` - original training data (currently unused beyond being
passed through); `validation_depths` - depths at which to test
predictions. **Returns**: A list with `model_validation` (per-property,
per-group prediction/diagnostic results), `prediction_validation`
(unused/reserved), `overall_assessment` (`total_models`, `valid_models`,
`success_rate`, `validation_passed` = success rate \>= 0.8), and
`recommendations`. **Algorithm**: For every property/group model,
predicts at `validation_depths` via
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
then checks `assess_trend_monotonicity()` (\>=70% of depth-to-depth
diffs share the same sign) and `assess_realistic_values()` (predictions
inside a property-specific realistic range table, e.g. clay/sand/silt in
`[0,100]`, pH in `[2.5,11]`, bulk density in `[0.3,3.0]`), logging
recommendations when models fail or produce unrealistic values.

**Minor exported/internal helpers** (not separately detailed above):
`create_soil_groups()` (assigns `soil_group` from the chosen taxonomy
column, `"unknown"` if missing),
`create_test_grouping()`/`apply_hierarchical_grouping()`/`apply_matching_hierarchy()`
(grouping and matching mechanics used by the functions above),
`generate_gp_processing_summary()` (logging only),
`create_model_summary()` (builds the `model_summary` metadata block),
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
and
[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)
(k-fold CV across three GPfit correlation families - exponential power
1.95, Matern nu=3/2, Matern nu=5/2 - picking the lowest mean held-out
RMSE; CV results are attached to the winning model as a `cv_results`
attribute), `calculate_model_diagnostics()` (training RMSE and
log-likelihood),
`assess_trend_monotonicity()`/`assess_realistic_values()` (used by
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)),
[`infer_simulation_properties()`](https://jjmaynard.github.io/soilSIM/reference/infer_simulation_properties.md)
(numeric, non-structural columns of a long-format simulation data
frame),
[`apply_nrcs_gp_adjustments_with_correlations()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_gp_adjustments_with_correlations.md)
/
[`apply_local_gp_adjustments_with_correlations()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments_with_correlations.md)
(thin wrappers that infer the `properties` argument and delegate to
`multivariate-adjustment.R`’s
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md)
/
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)),
`validate_enhanced_simulations()` (row-count-match check used by
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)),
[`select_simulation_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/select_simulation_correlation_matrix.md)
and
[`flatten_simulation_array_to_long()`](https://jjmaynard.github.io/soilSIM/reference/flatten_simulation_array_to_long.md)
(Monte Carlo input-shaping helpers for
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)),
and `export_gp_models()`/`load_gp_models()` (RDS-based persistence of a
fitted GP model set, optionally stripping training data to shrink file
size).

### multivariate-adjustment.R

#### `integrate_monte_carlo_with_gp()` - Master Function

``` r

integrate_monte_carlo_with_gp(
  simulation_results,
  gp_models = NULL,
  cokey_mapping = NULL,
  integration_method = "hybrid",
  preserve_correlations = TRUE,
  properties = NULL,
  parallel = FALSE,
  n_cores = NULL,
  config = NULL
)
```

**Purpose**: Whole-dataset entry point that applies GP depth-trend
adjustment (NRCS-model-based and/or locally-fit) to every cokey in a
Monte Carlo simulation result, with validation and rich metadata.
**Parameters**: `simulation_results` - output of
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)
(must have a `simulation_data` element); `gp_models` - fitted NRCS GP
models from
[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md)
(optional); `cokey_mapping` - output of
[`match_soils_to_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/match_soils_to_gp_models.md)
(optional); `integration_method` - `"nrcs_gp"` (NRCS models only),
`"local_gp"` (per-cokey local GP only), or `"hybrid"` (NRCS adjustment
followed by local GP adjustment); `preserve_correlations` - whether to
run the correlation-preserving quantile transform; `properties` -
properties to adjust, or `NULL` to auto-detect via
`detect_simulation_properties()`; `parallel`/`n_cores` - parallel
dispatch across cokeys; `config` - Module 8-style config, defaults to
`get_default_configuration("full")`. **Returns**: A list:
`integrated_data` (final adjusted long-format data frame),
`original_simulation_data`, `integration_metadata` (method, properties
processed, success rate, processing time, package versions, etc.),
`validation_results` (from
[`validate_integration_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_integration_results.md)),
`original_metadata`; also carries `success_rate`, `processing_time`, and
`validation_passed` as attributes. **Algorithm**: Validates parameters
and simulation-data quality, auto-detects/validates properties, decides
`use_nrcs_gp`/`use_local_gp` from `integration_method` and whether GP
models/mapping were supplied, then dispatches per-cokey work to
`process_cokeys_parallel()` or `process_cokeys_sequential()` (each
calling `process_single_cokey()`, which applies
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md)
then/or
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
per cokey inside a `tryCatch`), combines successful per-cokey results
with `combine_and_validate_results()`, runs
[`validate_integration_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_integration_results.md),
and packages everything via `create_integration_results()`.

#### `apply_gp_depth_trends()`

``` r

apply_gp_depth_trends(
  cokey_data,
  gp_predictions,
  properties,
  preserve_correlations = TRUE,
  primary_property = NULL
)
```

**Purpose**: Core single-cokey routine that reshapes long-format
simulation rows into matrices, applies the GP-trend adjustment
(correlation-preserving or independent), and reshapes back.
**Parameters**: `cokey_data` - long-format rows for one cokey;
`gp_predictions` - named list of depth-indexed GP mean predictions per
property; `properties` - properties to adjust; `preserve_correlations` -
use
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)
(\>=2 properties) vs. `apply_individual_adjustments()`;
`primary_property` - reference property for quantile ordering (defaults
to the first available property). **Returns**: `cokey_data` with
adjusted property columns merged back in (unchanged if there are fewer
than 2 rows, fewer than 2 valid depths, or no properties overlap with
`gp_predictions`). **Algorithm**: Filters to valid depths, builds
`[depth x simulation]` matrices via
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md),
adjusts them via
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)
or `apply_individual_adjustments()`, converts the result back to long
format via
[`convert_to_long_format()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_long_format.md),
and merges the adjusted values into `cokey_data` via
`merge_adjusted_data()` (matched on `hzdept_r` + `simulation_number`).
Every stage is wrapped in
`tryCatch`/[`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md)
so a failure at any step degrades to returning the
unmodified/partially-modified data rather than aborting.

#### `preserve_correlation_structure()`

``` r

preserve_correlation_structure(property_matrices, gp_predictions, depths, primary_property)
```

**Purpose**: The production correlation-preserving depth-adjustment
algorithm (long-format pipeline’s counterpart to `gp-modeling.R`’s
[`adjust_multivariate_depthwise_GP()`](https://jjmaynard.github.io/soilSIM/reference/adjust_multivariate_depthwise_GP.md)).
**Parameters**: `property_matrices` - named list of
`[depth x simulation]` matrices; `gp_predictions` - named list of
per-depth GP mean vectors; `depths` - depth vector; `primary_property` -
reference property (falls back to the first property if missing from
`property_matrices`). **Returns**: Named list of adjusted matrices, same
shape as `property_matrices`. **Algorithm**: Computes an ECDF of the
primary property’s surface-depth (`depths[1]`) values to fix
per-simulation `reference_quantiles`; for each property with a matching,
depth-length GP prediction vector, walks depths 2..n computing a
safety-clamped GP ratio (`calculate_safe_gp_ratio()`, clamped to
`[0.1, 10]`), nudges each simulation toward
`previous_adjusted_value * gp_ratio` at its reference quantile within
the current property’s own current-depth distribution
(`apply_quantile_adjustment()`), then remaps the nudged values back onto
the current property’s original ECDF (`correct_distribution_shape()`) to
preserve its marginal shape. Properties lacking a usable GP prediction
are passed through unchanged.

#### `apply_nrcs_trend_adjustments()`

``` r

apply_nrcs_trend_adjustments(cokey_data, gp_models, model_group, properties, preserve_correlations = TRUE)
```

**Purpose**: Adjust one cokey’s simulated properties toward
NRCS-model-derived depth trends for its matched GP model group.
**Parameters**: `cokey_data` - simulation rows for one cokey;
`gp_models` - fitted NRCS GP models; `model_group` - the GP model group
name for this cokey (from
[`match_soils_to_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/match_soils_to_gp_models.md)
/
[`match_simulations_to_nrcs_models()`](https://jjmaynard.github.io/soilSIM/reference/match_simulations_to_nrcs_models.md));
`properties` - simulation property names to adjust;
`preserve_correlations` - passed through to
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md).
**Returns**: Adjusted `cokey_data` (unchanged if no properties map to an
NRCS property name or fewer than 2 unique depths exist). **Algorithm**:
Maps simulation property names to NRCS property names via
`get_nrcs_property_mapping()` (e.g. `"sandtotal"` -\> `"sand_pct"`,
`"ph"` -\> `"pH"`, `"om"`/`"cec"`/`"soc"` -\> `"organic_matter"`),
obtains per-depth GP predictions for the matched group via
`get_nrcs_gp_predictions()` (falling back to local per-depth means when
no matching NRCS model exists for a property/group), then calls
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md).

#### `match_simulations_to_nrcs_models()`

``` r

match_simulations_to_nrcs_models(cokey, cokey_mapping, fallback_group = "general_pool")
```

**Purpose**: Look up a single cokey’s NRCS GP model group from a
precomputed mapping. **Returns**: The matched `gp_model_group` string,
or `fallback_group` if `cokey_mapping` is `NULL`, no row matches, or the
match is `NA`.

#### `extract_nrcs_depth_trends()`

``` r

extract_nrcs_depth_trends(gp_models, properties, depths = seq(0, 200, by = 10))
```

**Purpose**: Produce tidy depth-trend curves for every (property, group)
model, for reporting/plotting. **Returns**: A named list (by property)
of data frames with columns `depth`, `predicted_value`, `group`,
`property`, row-bound across all groups via
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html).
**Algorithm**: For each `"stratified_grouped"` property, calls
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md)
per group at the requested `depths` and assembles the results.

#### `apply_local_gp_adjustments()`

``` r

apply_local_gp_adjustments(
  cokey_data,
  properties,
  preserve_correlations = TRUE,
  min_depths = 3,
  config = NULL
)
```

**Purpose**: Adjust one cokey’s simulated properties toward depth trends
fitted *locally* from that cokey’s own data (used when no NRCS GP
model/mapping is available, or in `"local_gp"`/`"hybrid"` integration
mode). **Parameters**: `min_depths` - minimum unique depths required to
attempt local GP fitting (default 3); others as above. **Returns**:
Adjusted `cokey_data` (unchanged if too few unique depths, no local
models fit, or no valid predictions generated). **Algorithm**: Fits
per-property local GP models via
[`fit_local_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_models.md),
generates predictions at the cokey’s unique depths via
`generate_local_predictions()`, and applies them via
[`apply_local_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_depth_trends.md)
(which itself calls
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)).

#### `fit_local_gp_models()`

``` r

fit_local_gp_models(cokey_data, properties, config = NULL)
```

**Purpose**: Fit a small 1-D GP per property using only the current
cokey’s own depth-aggregated data. **Returns**: Named list of fitted
model objects (same shape as
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)’s
`model` element), only for properties with \>= 3 aggregated depth points
and non-zero variance. **Algorithm**: For each property, aggregates via
`aggregate_property_by_depth()`, checks variance, and fits via
`fit_local_gp_model_single()` (min-max depth scaling +
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)).

#### `apply_local_depth_trends()`

``` r

apply_local_depth_trends(cokey_data, local_predictions, unique_depths, preserve_correlations = TRUE)
```

**Purpose**: Thin wrapper applying locally-fit GP predictions via the
same
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)
machinery used for NRCS predictions.

#### `convert_to_property_matrices()`

``` r

convert_to_property_matrices(simulation_data, properties, unique_depths, sim_numbers)
```

**Purpose**: Reshape a long-format simulation data frame into one
`[depth x simulation]` matrix per property. **Returns**: Named list of
matrices (only for properties with at least one non-`NA` value found).
**Algorithm**: Nested loop over `unique_depths` x `sim_numbers`, pulling
the matching value via
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)/[`dplyr::pull()`](https://dplyr.tidyverse.org/reference/pull.html)
for each cell (documented as Module 8 “safe” lookups; O(depths x sims)
filters per property).

#### `convert_to_long_format()`

``` r

convert_to_long_format(adjusted_matrices, unique_depths, sim_numbers, original_data, properties)
```

**Purpose**: Inverse of
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md) -
reshape adjusted matrices back into a long-format data frame carrying
original row metadata (`cokey`, `compname`, `mukey`, `hzdept_r`,
`hzdepb_r`, `simulation_number`, `unique_id`). **Returns**: Long-format
data frame with one row per depth x simulation combination that had a
matching original row.

#### `validate_integration_results()`

``` r

validate_integration_results(original_data, integrated_data, properties, preserve_correlations = TRUE)
```

**Purpose**: Post-hoc quality check comparing pre- and post-integration
data. **Returns**: A list with `data_integrity` (row-count preservation,
missing-value increase per property, from `validate_data_integrity()`),
`correlation_preservation` (max/mean correlation differences at up to 3
depths, from `validate_correlation_preservation_integration()`, only
computed when `preserve_correlations` and \>=2 properties),
`trend_realism` (extreme-jump and out-of-realistic-range checks per
property, from `validate_depth_trends()`), and `overall_assessment`
(`integrity_score`, `correlation_score`, `trend_score`, `overall_score`
= their mean, `validation_passed` = `overall_score >= 0.8`).

#### `correct_distribution_shapes()`

``` r

correct_distribution_shapes(adjusted_data, original_data, properties, config = NULL)
```

**Purpose**: Post-adjustment cleanup pass that clamps each property to a
realistic range, optionally remaps its distribution shape back toward
the original, and enforces cross-property constraints (texture
sum-to-100). **Algorithm**: For each property, looks up
`get_property_constraints()` (range + whether to preserve distribution
shape), applies `apply_range_constraints()` (clamp), applies
`correct_property_distribution()` (ECDF quantile remap back onto
`original_data`) when the constraint calls for it, then applies
`apply_cross_property_constraints()` (rescales sand/silt/clay-family
columns so each row’s sum is exactly 100 whenever their sum is \> 0).

**Minor exported/internal helpers** (not separately detailed above):
`detect_simulation_properties()` (intersects
`get_predefined_properties("laboratory")` and a hardcoded
simulation-column-name list with the data’s actual columns),
`process_cokeys_parallel()` / `process_cokeys_sequential()` /
`process_single_cokey()` (cokey-level dispatch used by
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md) -
the parallel path uses
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html) +
`clusterEvalQ(cl, library(soilSIM))` on Windows so worker processes have
`dplyr`/`GPfit` loaded, or
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
elsewhere), `combine_and_validate_results()` /
`create_integration_results()` (result assembly and metadata),
`get_nrcs_property_mapping()` (the simulation-name -\> NRCS-name lookup
table), `get_nrcs_gp_predictions()` / `get_local_property_means()`
(per-property prediction retrieval with fallback),
`calculate_safe_gp_ratio()` (clamped consecutive-depth GP ratio),
`apply_quantile_adjustment()` / `correct_distribution_shape()` (the
per-depth quantile-nudge and ECDF-remap steps inside
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)),
`generate_local_predictions()`
([`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md)
over a set of locally-fit models), `aggregate_property_by_depth()` /
`fit_local_gp_model_single()` (used by
[`fit_local_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_models.md)),
`merge_adjusted_data()` (row-matched merge-back used by
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)),
`validate_data_integrity()` /
`validate_correlation_preservation_integration()` /
`validate_depth_trends()` /
`calculate_overall_validation_score_integration()` /
`log_validation_summary()` (the components of
[`validate_integration_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_integration_results.md)),
`get_property_constraints()` / `apply_range_constraints()` /
`correct_property_distribution()` / `apply_cross_property_constraints()`
(components of
[`correct_distribution_shapes()`](https://jjmaynard.github.io/soilSIM/reference/correct_distribution_shapes.md)),
`check_missing_values_increase()`, and `apply_individual_adjustments()`
(the non-correlation-preserving fallback path used by
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)
when `preserve_correlations = FALSE` or fewer than 2 properties are
present - scales each property independently by the clamped GP ratio
with no cross-property quantile linkage).

## Internal Connections

    build_stratified_gp_models() [MASTER, gp-modeling.R]
    ├── apply_hierarchical_grouping()
    │   └── create_test_grouping() [via select_optimal_grouping(), if strategy = "auto"]
    ├── fit_individual_gp_model()  (per property x final_group)
    │   ├── optimize_gp_hyperparameters()
    │   │   └── k_fold_gp_cv() -> GPfit::GP_fit() / GPfit::predict.GP()
    │   ├── GPfit::GP_fit()        [if optimize_hyperparameters = FALSE]
    │   └── calculate_model_diagnostics() -> GPfit::predict.GP()
    └── create_model_summary()

    prepare_nrcs_training_data() [gp-modeling.R]
    ├── is_unsuitable(), standardize_property_names(), safe_coalesce()  [external utility module]
    ├── select_optimal_grouping()  [if grouping_strategy = "auto"]
    ├── create_soil_groups()
    └── generate_gp_processing_summary()

    predict_gp_depth_trends() [gp-modeling.R]  <- called from BOTH files:
    ├── validate_gp_models()                     [gp-modeling.R]
    ├── extract_nrcs_depth_trends()               [multivariate-adjustment.R]
    ├── get_nrcs_gp_predictions()                 [multivariate-adjustment.R]
    └── generate_local_predictions()              [multivariate-adjustment.R]

    match_soils_to_gp_models() [gp-modeling.R]
    └── apply_matching_hierarchy()
            (consumed downstream by match_simulations_to_nrcs_models() in multivariate-adjustment.R)

    simulate_soil_properties() [gp-modeling.R, single-cokey convenience path]
    ├── select_simulation_correlation_matrix()
    ├── generate_monte_carlo_realizations()        [Monte Carlo module]
    ├── flatten_simulation_array_to_long()
    ├── apply_nrcs_gp_adjustments_with_correlations()
    │   └── apply_nrcs_trend_adjustments()         [-> multivariate-adjustment.R]
    ├── apply_local_gp_adjustments_with_correlations()
    │   └── apply_local_gp_adjustments()           [-> multivariate-adjustment.R]
    └── validate_enhanced_simulations()

    integrate_monte_carlo_with_gp() [MASTER, multivariate-adjustment.R, whole-dataset path]
    ├── detect_simulation_properties()
    ├── process_cokeys_parallel() / process_cokeys_sequential()
    │   └── process_single_cokey()                 (per cokey)
    │       ├── match_simulations_to_nrcs_models()
    │       ├── apply_nrcs_trend_adjustments()
    │       │   ├── get_nrcs_property_mapping()
    │       │   ├── get_nrcs_gp_predictions()
    │       │   │   └── predict_gp_depth_trends()  [-> gp-modeling.R]
    │       │   └── apply_gp_depth_trends()
    │       │       ├── convert_to_property_matrices()
    │       │       ├── preserve_correlation_structure() / apply_individual_adjustments()
    │       │       │   ├── calculate_safe_gp_ratio()
    │       │       │   ├── apply_quantile_adjustment()
    │       │       │   └── correct_distribution_shape()
    │       │       ├── convert_to_long_format()
    │       │       └── merge_adjusted_data()
    │       └── apply_local_gp_adjustments()
    │           ├── fit_local_gp_models() -> fit_local_gp_model_single() -> GPfit::GP_fit()
    │           ├── generate_local_predictions() -> predict_gp_depth_trends() [-> gp-modeling.R]
    │           └── apply_local_depth_trends() -> apply_gp_depth_trends()
    ├── combine_and_validate_results()
    ├── validate_integration_results()
    │   ├── validate_data_integrity()
    │   ├── validate_correlation_preservation_integration()
    │   ├── validate_depth_trends() -> get_property_constraints()
    │   └── calculate_overall_validation_score_integration()
    └── create_integration_results()

### Full GP + Integration Pipeline

    NRCS Combined Data (SSURGO components/horizons)
    ↓
    prepare_nrcs_training_data() [gp-modeling.R]
    ├── quality validation, name standardization, unsuitable-horizon filtering
    ├── auto grouping strategy selection [select_optimal_grouping()]
    └── property coalescing + soil_group assignment
    ↓
    build_stratified_gp_models() [gp-modeling.R]
    ├── hierarchical fallback grouping [apply_hierarchical_grouping()]
    └── per-group GP fitting [fit_individual_gp_model() -> GPfit::GP_fit()]
    ↓
    GP Models (by property x soil group) ──────────────┐
                                                         │
    Monte Carlo Simulation Output                       │
    (generate_monte_carlo_realizations(), separate module)
    ↓                                                    │
    match_soils_to_gp_models() [gp-modeling.R] ◄─────────┘
    ↓ (cokey -> gp_model_group mapping)
    integrate_monte_carlo_with_gp() [multivariate-adjustment.R, MASTER]
    ├── per-cokey NRCS GP adjustment  [apply_nrcs_trend_adjustments()]
    ├── per-cokey local GP adjustment [apply_local_gp_adjustments()]
    │   └── both funnel through apply_gp_depth_trends()
    │       ├── reshape to matrices [convert_to_property_matrices()]
    │       ├── correlation-preserving quantile adjustment
    │       │   [preserve_correlation_structure()]
    │       └── reshape back to long format [convert_to_long_format()]
    ↓
    combine_and_validate_results() + validate_integration_results()
    ↓
    Final Results: GP-depth-adjusted, correlation-preserving realizations per cokey
    ↓
    Statistics & Diagnostics (validation, summarization) [downstream module]

## Key Integration Points

### 1. Correlation-preserving quantile nudge (the core algorithm, both files)

``` r

# preserve_correlation_structure() / adjust_multivariate_depthwise_GP()
gp_ratio <- calculate_safe_gp_ratio(gp_means, i)          # clamped to [0.1, 10]
adjusted_curr <- apply_quantile_adjustment(
  reference_quantiles, curr_values, prev_values, gp_ratio, n_sims
)
adjusted_matrix[i, ] <- correct_distribution_shape(curr_values, adjusted_curr)
```

### 2. GP model group lookup feeding NRCS adjustment

``` r

# process_single_cokey() [multivariate-adjustment.R]
model_group <- match_simulations_to_nrcs_models(cokey, cokey_mapping)
result_data <- apply_nrcs_trend_adjustments(
  result_data, gp_models, model_group, properties, preserve_correlations
)
```

### 3. Cross-file delegation (gp-modeling.R wrappers -\> multivariate-adjustment.R)

``` r

# apply_nrcs_gp_adjustments_with_correlations() [gp-modeling.R]
apply_nrcs_trend_adjustments(
  cokey_data = simulation_data,
  gp_models = nrcs_gp_models,
  model_group = model_group,
  properties = infer_simulation_properties(simulation_data),
  preserve_correlations = preserve_correlations
)
```

### 4. Depth-scaled GP prediction shared by both files

``` r

# predict_gp_depth_trends() [gp-modeling.R] - called from get_nrcs_gp_predictions(),
# generate_local_predictions(), and extract_nrcs_depth_trends() in multivariate-adjustment.R
scaled_depths <- (new_depths - scaling$min) / scaling$range
predictions <- GPfit::predict.GP(gp_model_info$gp_model, xnew = as.matrix(scaled_depths))
```

## Dependencies

**External packages**: - `GPfit` -
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) /
[`GPfit::predict.GP()`](https://rdrr.io/pkg/GPfit/man/predict.html) are
the actual GP regression engine used throughout both files (single
covariate: scaled depth). - `parallel` -
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)/`clusterExport()`/`clusterEvalQ()`/`parLapply()`
on Windows,
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
elsewhere, used by `process_cokeys_parallel()` for cokey-level parallel
dispatch in
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md).
No `future`/`furrr` usage was found in either file - parallelism is base
`parallel` only. - `dplyr` - pervasive for
grouping/summarizing/filtering/joining throughout both files. - Base
`stats` - [`ecdf()`](https://rdrr.io/r/stats/ecdf.html),
[`quantile()`](https://rdrr.io/r/stats/quantile.html),
[`cor()`](https://rdrr.io/r/stats/cor.html),
[`var()`](https://rdrr.io/r/stats/cor.html),
[`sd()`](https://rdrr.io/r/stats/sd.html),
[`approx()`](https://rdrr.io/r/stats/approxfun.html) implement the
quantile-nudge/distribution-correction machinery.

**soilSIM dependencies (inputs)**: - Monte Carlo simulation output
([`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md),
its `simulation_data` array/long-format result) is the primary upstream
input to both
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
(single-cokey) and
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)
(whole-dataset). - Shared/general utility functions used throughout
(validation, logging, configuration, error handling) -
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md),
[`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md),
[`validate_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_parameters.md),
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
[`validate_properties()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties.md),
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md),
[`standardize_property_names()`](https://jjmaynard.github.io/soilSIM/reference/standardize_property_names.md),
[`safe_coalesce()`](https://jjmaynard.github.io/soilSIM/reference/safe_coalesce.md),
[`get_default_configuration()`](https://jjmaynard.github.io/soilSIM/reference/get_default_configuration.md),
[`get_predefined_properties()`](https://jjmaynard.github.io/soilSIM/reference/get_predefined_properties.md),
[`track_progress()`](https://jjmaynard.github.io/soilSIM/reference/track_progress.md),
`get_package_versions()` - these are not defined in either file and are
assumed to come from a shared validation/utility module elsewhere in the
package. - No direct call into a `distributions.R`-style
correlation-matrix-generation utility was found in either file;
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)/[`adjust_multivariate_depthwise_GP()`](https://jjmaynard.github.io/soilSIM/reference/adjust_multivariate_depthwise_GP.md)
preserve correlation implicitly via the shared reference-quantile
mechanism rather than by explicitly re-applying a target correlation
matrix (e.g. no Cholesky step appears in these two files - that logic
lives upstream, in the Monte Carlo generation step that
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
calls into via
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)).

**Downstream consumers**: - The `integrated_data`/`validation_results`
produced by
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)
feed the Statistics & Diagnostics group (percentile/summary statistics
computation and the broader validation-diagnostics workflow), which
consumes the GP-adjusted realizations as its input dataset and the
`validation_results` structure as a data-quality baseline.

## Data Flow In/Out

**In**: - NRCS/SSURGO combined component-horizon data (`cokey`,
`hzdept_r`/`hzdepb_r`, taxonomic columns, raw property columns) -\>
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md). -
Monte Carlo realizations (long-format simulation data frame with
`cokey`, `hzdept_r`, `simulation_number`, and per-property columns; or
the raw `[horizon, property, realization]` array from
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md))
-\>
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
/
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md). -
Fitted GP models
([`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md)
output) and a cokey-to-group mapping
([`match_soils_to_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/match_soils_to_gp_models.md)
output) -\> both integration entry points. - Correlation matrices (flat
or genhz-keyed) -\>
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
(consumed only to pick which matrix to hand to the Monte Carlo
generator; not directly re-applied in this file).

**Out**: - A fitted GP model set
([`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md)’s
return value: per-property, per-group models + diagnostics + summary),
typically persisted via `export_gp_models()`/`load_gp_models()`. -
GP-depth-adjusted, correlation-preserving realizations per cokey -
either as a single-cokey long data frame with validation attributes
([`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md))
or a whole-dataset result list
([`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)’s
`integrated_data` + `integration_metadata` + `validation_results`),
ready for percentile/statistical summarization downstream.

## Usage Example

``` r

# 1. Prepare NRCS training data and fit stratified GP depth-trend models
nrcs_training <- prepare_nrcs_training_data(
  nrcs_combined_data,
  grouping_strategy = "auto",
  max_depth = 200
)

gp_models <- build_stratified_gp_models(
  nrcs_training,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  optimize_hyperparameters = TRUE
)

# 2. Map the Monte Carlo simulation's cokeys onto GP model groups
cokey_mapping <- match_soils_to_gp_models(
  simulated_cokeys = unique(mc_results$simulation_data$cokey),
  nrcs_combined_data = nrcs_combined_data,
  gp_models = gp_models,
  property = "clay_pct"
)

# 3. Integrate the Monte Carlo output with the fitted GP models, preserving correlations
integrated <- integrate_monte_carlo_with_gp(
  simulation_results = mc_results,
  gp_models = gp_models,
  cokey_mapping = cokey_mapping,
  integration_method = "hybrid",
  preserve_correlations = TRUE,
  parallel = TRUE
)

final_realizations <- integrated$integrated_data
```
