# Monte Carlo Simulation

## Overview
`R/monte-carlo.R` is the core correlated Monte Carlo simulation engine of `soilSIM`. It takes SSURGO (or other) soil property data expressed as low/representative/high (`_l`/`_r`/`_h`) percentile triplets - already fitted to per-property distribution families by the Distribution Fitting core group in `distributions.R` - and produces many correlated realizations of soil properties per horizon. It handles both ordinary independent/correlated numeric properties (e.g. `dbovendry`, `ph1to1h2o`) and the compositional texture group (`sandtotal`/`silttotal`/`claytotal`), which is simulated jointly through an isometric-log-ratio (ILR) reparameterization so the sum-to-100% constraint holds exactly by construction rather than by post-hoc rescaling. The engine also optionally fuses caller-supplied observed field/lab data into each property's SSURGO-derived prior via Bayesian updating before simulating, applies a battery of physical/relationship/sum constraints after simulating, and returns rich validation, diagnostics, and quality-assessment metadata alongside the simulated values.

## Core Functions

Functions are grouped in the pipeline order that `generate_monte_carlo_realizations()` (the master orchestrator) calls them.

### 1. `validate_monte_carlo_config()` - Configuration Validation
**Purpose**: Validates a merged Monte Carlo configuration against a parameter spec before any data processing begins.
**Signature**:
```r
validate_monte_carlo_config(config, n_realizations)
```
**Parameters**:
- `config` - Configuration list to validate (nested under `$monte_carlo` or flat; both are handled).
- `n_realizations` - Number of realizations requested (used for an additional large-N/no-parallel warning).

**Returns**: List from `validate_parameters()` (`valid`, `errors`, `warnings`, plus per-field detail), with an extra warning appended if `n_realizations > 50000` and `parallel` is not enabled.
**Behavior**: Extracts `config$monte_carlo` if nested, then checks `distribution_type` (one of `triangular`/`normal`/`uniform`/`beta`/`lognormal`/`metalog`/`linear_cdf`/`auto`), `max_depth` (10-1000), `min_quality_score`/`min_success_rate`/`max_outlier_rate` (0-1), `outlier_threshold` (1-5), `error_recovery_action` (`stop`/`warn`/`continue`/`retry`), `normalization_method` (`proportional`/`additive`/`multiplicative`), and `correlation_fallback` (`identity`/`kssl_global`) via `validate_parameters()`.

### 2. `get_monte_carlo_defaults()` - Default Configuration
**Purpose**: Returns the full default Monte Carlo configuration, merged (safely) on top of the package's general default configuration.
**Signature**:
```r
get_monte_carlo_defaults()
```
**Returns**: A list with a `$monte_carlo` element containing every tunable default, including:
- `distribution_type = "triangular"` - global fallback family when a property has no `property_distributions` entry.
- `lh_percentile = c(0.05, 0.95)` - the assumed percentile meaning of SSURGO `_l`/`_h`.
- `property_distributions = list()` - per-property family/bounds overrides, e.g. `list(claytotal = list(family = "beta", bounds = c(0, 100)))`.
- `composition_groups = list(texture = list(members = c("sandtotal", "silttotal", "claytotal"), pseudo = c("ilr1", "ilr2")))` - the sand/silt/clay role order is chosen to match the sequential binary partition used to build the KSSL reference correlation matrices, so the `kssl_global` fallback's ilr1/ilr2 correlations are directly reusable (this reparameterization is a pure internal detail; simulated clay/sand/silt output is invariant to role order).
- `max_depth = 250`, `auto_correlation = FALSE`, `correlation_fallback = "identity"` (opt-in `"kssl_global"` alternative), `parallel_threshold = 1000`.
- Quality-control knobs: `min_quality_score = 0.6`, `min_success_rate = 0.8`, `outlier_threshold = 2.5`, `max_outlier_rate = 0.1`, `min_data_points = 5`, `minimum_observations = 10`.
- Error handling: `error_recovery_action = "warn"`, `strict_property_validation = FALSE`, `strict_unsuitability = FALSE`.
- Constraint toggles: `apply_range_constraints`/`apply_sum_constraints`/`apply_relationship_constraints = TRUE`, `apply_physical_constraints = FALSE`.
- Component simulation: `normalization_method = "proportional"`, `composition_constraints = NULL`.
- `quality_thresholds` list (`completeness_threshold`, `missing_data_threshold`, `outlier_threshold`, `correlation_threshold`).

### 3. `normalize_monte_carlo_config()` - Flat/Nested Config Normalization
**Purpose**: Accepts either a flat `simulation_config` (e.g. `list(distribution_type = "normal")`, the shape shown in `generate_monte_carlo_realizations()`'s own `@examples`) or an already-nested one (`list(monte_carlo = list(...))`), and normalizes to the nested shape `merge_configurations()` expects.
**Signature**:
```r
normalize_monte_carlo_config(simulation_config, default_config)
```
**Parameters**:
- `simulation_config` - User-supplied config, flat or nested.
- `default_config` - Result of `get_monte_carlo_defaults()`, used only to recognize which flat key names belong under `monte_carlo`.

**Returns**: `simulation_config`, wrapped under `$monte_carlo` if it was flat and its names matched known `monte_carlo` keys; passed through unchanged otherwise.
**Behavior/rationale**: A flat config previously landed silently as an unused top-level key (confirmed by testing: passing `distribution_type = "normal"` flat left `config$monte_carlo$distribution_type` at `"triangular"` every time) since `merge_configurations()` merges strictly by key path. This function is called immediately before `merge_configurations()` in `generate_monte_carlo_realizations()` to fix that.

### 4. `validate_monte_carlo_inputs()` - Input Validation
**Purpose**: Comprehensive validation of `soil_data`, `properties`, an optional `correlation_matrix`, and `n_realizations` before any simulation work happens.
**Signature**:
```r
validate_monte_carlo_inputs(soil_data, properties, correlation_matrix, n_realizations, config)
```
**Parameters**:
- `soil_data` - Input soil data frame (must carry `_l`/`_r`/`_h` columns per property, plus `cokey`).
- `properties` - Character vector of properties to simulate.
- `correlation_matrix` - Optional user-supplied correlation matrix.
- `n_realizations` - Requested realization count.
- `config` - Merged simulation configuration.

**Returns**: `list(valid=, errors=, warnings=, data_validation=, property_validation=, parameter_validation=)`.
**Behavior**: Runs `validate_parameters()` on `n_realizations` (range 1-100000); runs `validate_data_quality()` against required column `cokey` and the expected `_l/_r/_h` numeric columns, using safe (`tryCatch`-wrapped) access to `quality_thresholds` with hardcoded fallbacks; runs `validate_properties_with_synonyms()` against the SSURGO property lookup (warnings only unless `strict_property_validation`); validates any user `correlation_matrix` via `validate_correlation_matrix()` (warning only, not fatal); and calls `check_data_sufficiency()` to warn if too few properties have enough non-missing representative values.

### 5. `resolve_real_properties()` - Composition-Group Column Resolution
**Purpose**: Small internal helper resolving composition-group pseudo-properties (e.g. `"ilr1"`/`"ilr2"`) back to their real backing SSURGO column names (e.g. `sandtotal`/`silttotal`/`claytotal`) for data-availability checks, since a literal `"ilr1_r"` column never exists.
**Signature**:
```r
resolve_real_properties(properties, composition_plan = NULL)
```
**Parameters**: `properties` - possibly including pseudo-properties; `composition_plan` - result of `distributions.R`'s `resolve_composition_groups()`, or `NULL` (pass-through).
**Returns**: Character vector of real property names, deduplicated.

### 6. `prepare_simulation_data()` - Data Preparation
**Purpose**: Filters the input `soil_data` down to horizons suitable for simulation.
**Signature**:
```r
prepare_simulation_data(soil_data, properties, config, composition_plan = NULL)
```
**Parameters**:
- `soil_data` - Input soil data.
- `properties` - Properties to prepare, possibly including composition-group pseudo-properties.
- `config` - Simulation configuration.
- `composition_plan` - Optional result of `resolve_composition_groups()`; used via `resolve_real_properties()` so data-availability filtering checks the real underlying columns.

**Returns**: The filtered `soil_data` (only rows passing all masks).
**Behavior**: Computes `unsuitable_horizon` via `is_unsuitable()` if not already present; combines that with a depth filter (`hzdepb_r <= config$monte_carlo$max_depth`, when `hzdepb_r` exists) and a per-row data-availability mask from `check_property_data_availability()` (a row passes if it has non-missing `_r` data for at least one requested real property). Logs the row count before/after filtering.

### 7. `prepare_simulation_parameters()` - Parameter Extraction
**Purpose**: Builds the per-horizon, per-property fitted-distribution parameter list that drives simulation, special-casing active composition groups.
**Signature**:
```r
prepare_simulation_parameters(simulation_data, properties, config, composition_plan = NULL)
```
**Parameters**:
- `simulation_data` - Filtered simulation data (post `prepare_simulation_data()`).
- `properties` - Properties to extract parameters for (may include pseudo-properties).
- `config` - Simulation configuration.
- `composition_plan` - Optional result of `resolve_composition_groups()`.

**Returns**: A list of length `nrow(simulation_data)`, each element a named list keyed by property, each value `list(family=, fit=, source=)`.
**Behavior**: For "ordinary" (non-composition-group) properties, calls `extract_property_parameters()` per horizon per property. For each *active* composition group's pseudo-properties (e.g. `ilr1`/`ilr2`), both are fit together, once per horizon, via `distributions.R`'s `estimate_ilr_moments_mc()` from the group's real `_l/_r/_h` triplets (their joint covariance only makes sense computed jointly). Only the marginal SDs from that MC estimate are kept (`family = "normal"`, `source = "ilr_mc"`); the ilr1<->ilr2 *off-diagonal* is deliberately discarded because the later cross-horizon empirical correlation step (`estimate_property_correlations()`) owns that correlation instead - avoiding double-counting per-observation uncertainty propagation against cross-sample empirical correlation. Finishes by logging a per-property coverage summary via `summarize_parameter_extraction()`.

### 8. `extract_property_parameters()` - Single Horizon/Property Parameter Fit
**Purpose**: Resolves one property's distribution family and fits family-appropriate parameters from its SSURGO `_l/_r/_h` triplet for one horizon.
**Signature**:
```r
extract_property_parameters(horizon_data, property_name, config)
```
**Parameters**:
- `horizon_data` - One-row data frame/list for a single horizon.
- `property_name` - Property name without the `_l`/`_r`/`_h` suffix.
- `config` - Simulation configuration (consults `config$monte_carlo$property_distributions[[property_name]]` first, falling back to `config$monte_carlo$distribution_type`, default `"triangular"`).

**Returns**: `list(valid=, parameters=list(family=, fit=, source=))`.
**Behavior**: If all three of `l`/`r`/`h` are finite, delegates to `distributions.R`'s `fit_percentile_triplet()` with the resolved family, `lh_percentile`, and any registered `bounds` (`source = "complete"`). If only `r` is finite, falls back to a simple +/-20%-around-`r` triangular estimate, explicitly tagged `family = "triangular"` regardless of the requested family (`source = "estimated"`), so `quantile_from_fit()` still dispatches correctly. If `r` itself is missing, returns `valid = FALSE`. Prior to this being implemented for real, this function only ever returned `list(min=, mode=, max=)` regardless of requested family, so `"normal"`/`"beta"` requests silently produced downstream `NA`s.

### 9. `fuse_observed_data_into_priors()` - Bayesian Fusion of Observed Data (optional)
**Purpose**: Orchestrates Bayesian updating of each property's per-horizon prior (from `prepare_simulation_parameters()`) against a caller-supplied observed-data likelihood, using `bayesian-updating.R`'s pure fusion primitives.
**Signature**:
```r
fuse_observed_data_into_priors(simulation_params, simulation_data, sim_properties,
                                properties, observed_data, composition_plan, config)
```
**Parameters**:
- `simulation_params` - Per-horizon parameter list from `prepare_simulation_parameters()`.
- `simulation_data` - The filtered per-horizon data frame (same row order as `simulation_params`).
- `sim_properties` - Property vector `simulation_params` is keyed over (may include ILR pseudo-properties).
- `properties` - Caller-facing original property vector (kept for interface symmetry; not used directly).
- `observed_data` - Named list keyed by property name; texture entries are supplied jointly as `claytotal`/`sandtotal`/`silttotal`.
- `composition_plan` - Result of `resolve_composition_groups()`.
- `config` - Simulation configuration.

**Returns**: `simulation_params` with fused entries replaced in place (unfused entries pass through unchanged).
**Behavior**: Two parts. (1) **Texture (composition group) fusion**: if all three texture members are present in `observed_data`, each is reduced to a `c(low, rep, high)` triplet (via the internal `as_lrh_triplet()` - accepting named `low/rep/high` or `l/r/h` triplets, an unnamed length-3 vector, or a raw sample vector reduced via its own empirical quantiles at `lh_percentile`), and `fuse_texture_group_from_triplets()` is called once per horizon (that horizon's own clay/sand/silt `_l/_r/_h` as the prior side, the shared observed triplet as the likelihood side), replacing `simulation_params[[i]][["ilr1"]]`/`[["ilr2"]]`. Partial texture entries (1-2 of 3) skip texture fusion with a logged WARN. (2) **Ordinary property fusion**: for every non-texture property present in both `observed_data` and `sim_properties`, each horizon's prior is fused via the internal `fuse_one_property_prior()` helper against the *same* shared observed likelihood (every horizon keeps its own prior but shares one likelihood, since field/lab measurements aren't usually tied to individual simulated horizons). Two supported likelihood shapes: a raw numeric sample vector routes through the fully general `fuse_property()` grid-KDE path (prior is first sampled via `quantile_from_fit()`, posterior stored as an exact empirical quantile function with `family = "linear_cdf"`); a family-native parameter list (`list(mean=, sd=)` or `list(shape1=, shape2=)`/`list(mean=, var=)`) routes through the closed-form `bayes_fuse()`-based path, supported only for prior families `normal`/`lognormal`/`beta` (the only families with a conjugate route) - any other prior family skips fusion for that (horizon, property) with a logged WARN. Note: this only narrows each property's *marginal* prior spread; it does not touch the cross-property correlation/dependency structure, which `estimate_property_correlations()` estimates separately from raw representative values.

Minor related internal helpers (not exported, `@keywords internal`):
- `as_lrh_triplet(x, lh_probs)` - reduces an `observed_data` entry to an unnamed `c(low, rep, high)` triplet.
- `fuse_one_property_prior(prior, likelihood, is_vector_likelihood, n_samples, supported_closed_form, prop, horizon_index)` - implements the per-property closed-form/general dispatch described above, handling `normal`, `lognormal` (converts the raw-space likelihood to log-space via `normal_to_lognormal_params()` before fusing, since the lognormal prior fit is already stored in log-space), and `beta` (rescales a raw-scale likelihood mean/var onto the prior's own `[lower, upper]` support before fitting moments, since `fuse_beta()`'s alpha/beta addition requires shared support; falls back to a Normal-moments round-trip via `bayes_update_normal_normal()` if the direct beta fusion is infeasible).

### 10. `configure_correlation_structure()` - Correlation Matrix Setup
**Purpose**: Determines the correlation matrix to use for the Cholesky-copula simulation step, from a user-supplied matrix, auto-estimation, or independence.
**Signature**:
```r
configure_correlation_structure(simulation_params, properties, correlation_matrix = NULL,
                                 config, simulation_data = NULL)
```
**Parameters**:
- `simulation_params` - Simulation parameters (post fusion, if any).
- `properties` - Property names (drives matrix dimension/order).
- `correlation_matrix` - Optional user-supplied correlation matrix.
- `config` - Simulation configuration (`config$monte_carlo$auto_correlation`, `correlation_fallback`).
- `simulation_data` - Optional raw per-horizon data frame, passed through to `estimate_property_correlations()` so genhz-stratified estimation can be used when a `genhz` column is present.

**Returns**: `list(matrix=, method=, properties=, validation=, statistics=)`. `method` is one of `"user_provided"`, `"identity_fallback"`, `"estimated_from_data"`, `"identity_estimation_failed"`, or `"identity_default"`.
**Behavior**: If `correlation_matrix` is supplied, validates it via `validate_correlation_matrix()` and uses it if valid (falls back to identity with a WARN otherwise). Else if `auto_correlation` is `TRUE`, calls `estimate_property_correlations()`. Else uses the identity matrix (independent simulation). In all cases, sets row/column names to `properties`, repairs positive-definiteness via `ensure_positive_definite_matrix()`, computes matrix statistics (`calculate_matrix_statistics()`: eigenvalues, condition number, determinant), and does a final `validate_correlation_matrix()` pass.

### 11. `estimate_property_correlations()` - Empirical Correlation Estimation
**Purpose**: Builds a real empirical correlation matrix from the per-horizon fitted representative values (mean/mode per property, including composition pseudo-properties), optionally stratified by genetic horizon (`genhz`), with an opt-in KSSL reference-matrix fallback.
**Signature**:
```r
estimate_property_correlations(simulation_params, properties, config, simulation_data = NULL)
```
**Parameters**:
- `simulation_params` - Per-horizon parameter list.
- `properties` - Character vector of property names.
- `config` - Simulation configuration.
- `simulation_data` - Optional raw per-horizon data frame (same row order); its `genhz` column, if present, enables genhz-stratified estimation (explicit `genhz` always takes precedence over auto-derivation from `hzname`).

**Returns**: `list(valid=, correlation_matrix=, info=list(method=, n_obs=))`.
**Behavior**: Extracts one representative value per horizon per property via an internal `representative_value()` switch over distribution family (mode for triangular, midpoint for uniform, mean for normal, `exp(mean)` for lognormal, the Beta mean rescaled to `[lower, upper]`, the median for metalog/linear_cdf). Assembles these into a data frame and delegates to `distributions.R`'s `estimate_correlation_matrix_robust()` (`Hmisc::rcorr()` plus positive-definiteness repair), grouped by `genhz` when available. When `config$monte_carlo$correlation_fallback == "kssl_global"` (default `"identity"`), a group that fails empirical estimation falls back to `kssl-reference-correlations.R`'s `build_kssl_fallback_matrix()` for that group instead of being dropped, the global fallback becomes a KSSL-pooled matrix, and `genhz` is auto-derived from `hzname` via `classify_genhz()` if not already present. Reorders the final matrix to guarantee it exactly matches `properties`' order for the downstream Cholesky step. Prior to being implemented for real, this was a placeholder always returning `valid = FALSE`, silently forcing `auto_correlation = TRUE` back to the identity matrix.

### 12. `setup_distributions()` - Distribution Configuration Assembly
**Purpose**: Validates and packages each horizon/property's fitted parameters into a distribution-config object, with fallback handling for invalid fits.
**Signature**:
```r
setup_distributions(simulation_params, properties, config)
```
**Parameters**: `simulation_params` - per-horizon parameter list; `properties` - property names; `config` - simulation configuration (`config$monte_carlo$distribution_type` as the label fallback).
**Returns**: `list(distributions=, summary_stats=, validation_results=, configuration=)`, where `distributions` is a per-horizon list of per-property configs from `create_distribution_config()`/`create_fallback_distribution()`.
**Behavior**: For each horizon and property with existing parameters, calls `validate_distribution_parameters()` (a thin wrapper over `distributions.R`'s `validate_fit_parameters()`); on success wraps it via `create_distribution_config()`, on failure substitutes a `create_fallback_distribution()` (+/-20% triangular around the fit's mode). Then computes `summarize_distributions()` (counts/family breakdown) and `validate_distribution_setup()` (validity rate) over the whole set.

### 13. `simulate_correlated_properties()` - Core Simulation Kernel
**Purpose**: The actual correlated-random-number-generation engine: for every horizon, generates `n_realizations` correlated draws across `properties` using a Gaussian (Cholesky) copula and each property's fitted quantile function.
**Signature**:
```r
simulate_correlated_properties(simulation_params, correlation_matrix, n_realizations, properties, config)
```
**Parameters**:
- `simulation_params` - List containing simulation parameters for each horizon.
- `correlation_matrix` - Correlation matrix for the properties (identity used if `NULL`).
- `n_realizations` - Number of realizations to generate.
- `properties` - Character vector of property names.
- `config` - Configuration settings (used only incidentally here; each property's `family` travels with its own fitted `params`).

**Returns**: A `[horizon, property, realization]` numeric array (dimnamed accordingly).
**Behavior**: Initializes the result array; defaults/repairs the correlation matrix (`ensure_positive_definite_matrix()`, with a `validate_correlation_matrix()` check that falls back to identity on failure); Cholesky-decomposes it once (`chol()`). Per horizon: draws an independent uniform matrix (`n_properties x n_realizations`), transforms to standard normal (`qnorm`), applies the correlation structure (`chol_matrix %*% normal_matrix`), transforms back to correlated uniforms (`pnorm`), then for each property calls `distributions.R`'s `quantile_from_fit()` with that property's own resolved `family`/`fit` (family travels with params - `config$monte_carlo$distribution_type` is only ever consulted upstream as *their* fallback, never here) to get the final simulated values. Missing per-property parameters yield `NA` for that property/horizon. Reports progress via `track_progress()` and finishes with a `validate_simulated_array()` finite-rate check.

### 14. `apply_simulation_constraints()` - Post-Simulation Constraints
**Purpose**: Applies range, sum-to-100 (composition), relationship (ordering), and optional physical-plausibility constraints to the raw simulated array.
**Signature**:
```r
apply_simulation_constraints(simulation_results, properties, config, composition_plan = NULL)
```
**Parameters**:
- `simulation_results` - Array of simulation results (already composition-restored, i.e. no ILR pseudo-properties).
- `properties` - Property names.
- `config` - Configuration settings (`apply_range_constraints`/`apply_sum_constraints`/`apply_relationship_constraints`/`apply_physical_constraints` toggles).
- `composition_plan` - Optional result of `resolve_composition_groups()`; when the texture group is active, its sum-to-100 constraint is already exact by construction (via the ILR inverse in `restore_composition_properties()`), so `get_sum_constraints()` skips re-rescaling it (which would distort the correlation structure ILR exists to preserve).

**Returns**: The constrained array, with a `"constraint_summary"` attribute recording what was applied.
**Behavior**: Builds constraint rules via `get_constraint_rules()` (range: sand/silt/clay in `[0,100]`, `dbovendry` in `[0.3,3.0]`, pH in `[2.5,11.0]`; sum: any non-ILR texture subset re-rescaled to sum 100 with 5% tolerance; relationship: `wfifteenbar <= wthirdbar`, clipping violations; physical: optional wider soil-science bounds table for properties like `cec7`, `om`, `oc`, `rfv`, `awc`, etc., only for properties not already covered by range constraints), then sequentially applies `apply_range_constraints_batch()`, `apply_sum_constraints()`, `apply_relationship_constraints()`, and (if enabled) `apply_physical_constraints()`.

### 15. `validate_simulation_output()` - Output Validation
**Purpose**: Validates the final constrained simulation array's dimensions and per-property statistical quality.
**Signature**:
```r
validate_simulation_output(simulation_results, properties, config)
```
**Parameters**: `simulation_results` - the constrained array; `properties` - property names; `config` - configuration (`min_success_rate`, `outlier_threshold`, `max_outlier_rate`).
**Returns**: `list(valid=, success_rate=, property_validation=, quality_metrics=, issues=)`.
**Behavior**: First checks the array has 3 dimensions and the expected property count (fails fast otherwise). For each property, calls `analyze_property_simulation_quality()` (finite-rate, mean/SD, range) and flags a low success rate; runs `detect_outliers()` (IQR method) and flags a high outlier rate (guarding against `NaN` outlier rates - e.g. a property missing from every horizon - which previously crashed the pipeline via an `if (NaN > x)` error instead of degrading gracefully). Computes `calculate_simulation_quality_metrics()` for an overall quality summary, an overall `success_rate` across all properties/horizons/realizations, and sets `valid = FALSE` if that overall rate falls below `min_success_rate`.

### 16. `generate_monte_carlo_realizations()` - Master Pipeline
**Purpose**: The top-level, end-to-end orchestrator: validates inputs and config, prepares data and parameters, optionally fuses observed data, configures correlation and distributions, runs the simulation (sequential or parallel), restores compositional properties, applies constraints, validates output, and compiles diagnostics/quality assessment.
**Signature**:
```r
generate_monte_carlo_realizations(soil_data,
                                   properties,
                                   correlation_matrix = NULL,
                                   n_realizations = 1000,
                                   simulation_config = list(),
                                   validate_inputs = TRUE,
                                   parallel = FALSE,
                                   seed = NULL,
                                   observed_data = NULL)
```
**Parameters**:
- `soil_data` - Data frame with soil property data (must have `_l`/`_r`/`_h` columns per property, plus `cokey`).
- `properties` - Character vector of properties to simulate.
- `correlation_matrix` - Optional correlation matrix for the properties; if an explicit user-supplied matrix collides with an active composition group (a 3-variable raw sub-block has no principled remap onto a 2-variable ILR pair), the affected group falls back to independent per-property simulation with a logged WARN (auto-estimation has no such problem, since it computes its matrix from `sim_properties` directly).
- `n_realizations` - Integer, number of Monte Carlo realizations (default `1000`).
- `simulation_config` - List of simulation configuration parameters, flat or nested under `monte_carlo` (see `normalize_monte_carlo_config()`).
- `validate_inputs` - Logical, whether to run `validate_monte_carlo_inputs()` (default `TRUE`).
- `parallel` - Logical, whether to use parallel processing when `n_realizations` exceeds `config$monte_carlo$parallel_threshold` (default `FALSE`).
- `seed` - Integer, random seed for reproducibility (`set.seed()`).
- `observed_data` - Optional named list keyed by property name (plus, for the texture group, `claytotal`/`sandtotal`/`silttotal` supplied together) providing a Bayesian-updating likelihood to fuse against each property's SSURGO-derived prior before simulating - see `fuse_observed_data_into_priors()`. `NULL` (default) skips fusion entirely, preserving prior-only behavior exactly.

**Returns**: A list with `simulation_data` (final constrained `[horizon, property, realization]` array), `original_data` (the filtered pre-simulation data frame), `parameters` (per-horizon fitted parameters), `correlation_structure`, `distribution_setup`, `validation` (`input_validation`/`property_validation`/`output_validation`), `diagnostics`, `quality_assessment`, and `metadata` (properties, n_realizations, n_horizons, max_depth, processing_time, parallel, seed, timestamp, config_used, success_rate).
**Behavior**: See the numbered pipeline steps in "Internal Connections" below - this function is a straight-line orchestration of essentially every other function in this file, in a fixed step order (numbered 1 through 11 in the source's own logging), with two special orderings worth noting: (a) an explicit user `correlation_matrix` colliding with an active composition group causes `composition_plan`/`sim_properties`/`simulation_params` to be recomputed with that group forced inactive; (b) `fuse_observed_data_into_priors()` runs deliberately *after* that recompute (not immediately after parameter extraction), since the recompute replaces `simulation_params` wholesale and would otherwise silently discard fused posteriors.

### Other exported/minor helpers
- `sim_component_compositions(component_data, n_realizations = 1000, config = NULL)` - A separate, self-contained pipeline (not called by `generate_monte_carlo_realizations()`) for simulating SSURGO *component* percent composition (`comppct_l/_r/_h`), independent of the horizon-property engine above. Validates via `validate_component_data()`, extracts per-component min/mode/max via `extract_component_parameters()` (filling missing `comppct_l`/`comppct_h` from `comppct_r -/+ 2`), draws values via `generate_component_values()` (currently a uniform-distribution placeholder, not family-aware), normalizes each realization to sum to 100 via `normalize_component_realizations()`, optionally clamps via `apply_composition_constraints()`, and scores fidelity via `assess_component_quality()` (mean relative error of simulated vs. original `comppct_r`). Returns `list(realizations=, component_data=, n_realizations=, n_components=, quality_metrics=, validation=, constraints_applied=, metadata=)`.
- `run_sequential_simulation()` / `run_parallel_simulation()` - internal (not exported) dispatch helpers. The sequential path is a direct passthrough to `simulate_correlated_properties()`. The parallel path splits `n_realizations` into near-equal contiguous chunks across `parallel::detectCores() - 1` (or `config$monte_carlo$n_cores`) workers - `parallel::makeCluster()`/`clusterExport()`/`clusterEvalQ(library(soilSIM))`/`parLapply()` on Windows, `parallel::mclapply()` elsewhere - concatenates the resulting per-worker arrays along the realization dimension, and falls back to the sequential path on any error or if fewer than 2 usable cores are available.
- `calculate_matrix_statistics()`, `get_constraint_rules()`/`get_range_constraints()`/`get_sum_constraints()`/`get_relationship_constraints()`/`get_physical_constraints()`, `apply_range_constraints_batch()`/`apply_sum_constraints()`/`apply_relationship_constraints()`/`apply_physical_constraints()` - internal building blocks for `configure_correlation_structure()` and `apply_simulation_constraints()` respectively, described inline above.
- `check_data_sufficiency()` / `check_property_data_availability()` - internal per-property/per-row data-sufficiency checks feeding `validate_monte_carlo_inputs()` and `prepare_simulation_data()`.
- `analyze_property_simulation_quality()`, `calculate_simulation_quality_metrics()`, `generate_simulation_diagnostics()`, `assess_simulation_quality()` - internal post-simulation reporting helpers building `output_validation$property_validation`, `diagnostics`, and the final `quality_assessment` (an overall score blending `output_validation$success_rate`, mean-based "constraint satisfaction," and SD-based "statistical consistency").
- `validate_distribution_parameters()`, `create_distribution_config()`, `create_fallback_distribution()`, `summarize_distributions()`, `validate_distribution_setup()` - internal building blocks for `setup_distributions()`.
- `merge_configurations_safe()` - a `tryCatch`-guarded wrapper around `utils.R`'s `merge_configurations()` with a manual nested-list merge fallback if that function is somehow unavailable.

## Internal Connections

```
generate_monte_carlo_realizations()  [MASTER PIPELINE]
├── setup_logging()                                  [utils.R]        (if not already configured)
├── set.seed(seed)                                    (if seed supplied)
├── get_monte_carlo_defaults()
├── normalize_monte_carlo_config()
├── merge_configurations()                            [utils.R]
├── validate_monte_carlo_config()
│   └── validate_parameters()                         [utils.R]
├── Step 1: validate_monte_carlo_inputs()             (if validate_inputs)
│   ├── validate_parameters()                         [utils.R]
│   ├── validate_data_quality()                       [utils.R]
│   ├── validate_properties_with_synonyms()           [utils.R]
│   ├── validate_correlation_matrix()                 [distributions.R]
│   └── check_data_sufficiency()
├── Step 2: validate_properties_with_synonyms()       [utils.R]
├── Step 2.5: resolve_composition_groups()            [distributions.R]
├── Step 3: prepare_simulation_data()
│   ├── is_unsuitable()                               [utils.R]
│   ├── resolve_real_properties()
│   └── check_property_data_availability()
├── Step 4: prepare_simulation_parameters()
│   ├── extract_property_parameters()                 (ordinary properties)
│   │   └── fit_percentile_triplet()                  [distributions.R]
│   ├── estimate_ilr_moments_mc()                      [distributions.R]  (active composition groups)
│   └── summarize_parameter_extraction()
├── [conditional] recompute composition_plan/sim_properties/simulation_params
│   if an active group collides with a user-supplied correlation_matrix
├── Step 4.5: fuse_observed_data_into_priors()         (if observed_data supplied)
│   ├── as_lrh_triplet()
│   ├── fuse_texture_group_from_triplets()             [bayesian-updating.R]  (texture group)
│   └── fuse_one_property_prior()                      (ordinary properties)
│       ├── quantile_from_fit()                        [distributions.R]
│       ├── fuse_property()                            [bayesian-updating.R]
│       ├── bayes_update_normal_normal()                [bayesian-updating.R]
│       ├── normal_to_lognormal_params()                [bayesian-updating.R]
│       └── fuse_beta() / moments_to_beta() / beta_to_moments()  [bayesian-updating.R]
├── Step 5: configure_correlation_structure()
│   ├── validate_correlation_matrix()                  [distributions.R]  (user-provided path)
│   ├── estimate_property_correlations()                (auto_correlation path)
│   │   ├── estimate_correlation_matrix_robust()         [distributions.R]
│   │   ├── build_kssl_fallback_matrix()                 [kssl-reference-correlations.R]  (kssl_global fallback)
│   │   └── classify_genhz()                             [kssl-reference-correlations.R]  (kssl_global fallback)
│   ├── ensure_positive_definite_matrix()               [distributions.R]
│   └── calculate_matrix_statistics()
├── Step 6: setup_distributions()
│   ├── validate_distribution_parameters()
│   │   └── validate_fit_parameters()                   [distributions.R]
│   ├── create_distribution_config() / create_fallback_distribution()
│   ├── summarize_distributions()
│   └── validate_distribution_setup()
├── Step 7: run_sequential_simulation() or run_parallel_simulation()
│   └── simulate_correlated_properties()                [CORE KERNEL]
│       ├── ensure_positive_definite_matrix()           [distributions.R]
│       ├── validate_correlation_matrix()               [distributions.R]
│       ├── quantile_from_fit()                         [distributions.R]
│       └── validate_simulated_array()
├── Step 7.5: restore_composition_properties()          [distributions.R]  (ILR -> raw texture columns)
├── Step 8: apply_simulation_constraints()
│   ├── get_constraint_rules() -> get_range/_sum/_relationship/_physical_constraints()
│   └── apply_range/_sum/_relationship/_physical_constraints_batch()
├── Step 9: validate_simulation_output()
│   ├── analyze_property_simulation_quality()
│   ├── detect_outliers()                               [utils.R]
│   └── calculate_simulation_quality_metrics()
├── Step 10: generate_simulation_diagnostics()
└── Step 11: assess_simulation_quality()
```

## Dependencies

**From `distributions.R`** (percentile-triplet fitting, ILR compositional machinery, correlation-matrix utilities):
`fit_percentile_triplet()`, `quantile_from_fit()`, `validate_fit_parameters()`, `resolve_composition_groups()`, `restore_composition_properties()`, `estimate_ilr_moments_mc()`, `ensure_positive_definite_matrix()`, `validate_correlation_matrix()`, `estimate_correlation_matrix_robust()`.

**From `bayesian-updating.R`** (only exercised when `observed_data` is supplied to `generate_monte_carlo_realizations()`):
`bayes_fuse()`, `fuse_property()`, `fuse_texture_group_from_triplets()`, `bayes_update_normal_normal()`, `normal_to_lognormal_params()`, `fuse_beta()`, `moments_to_beta()`, `beta_to_moments()`.

**From `kssl-reference-correlations.R`** (only exercised when `config$monte_carlo$correlation_fallback == "kssl_global"`):
`build_kssl_fallback_matrix()`, `classify_genhz()`.

**From `utils.R`** (cross-cutting workflow infrastructure):
`log_message()`, `handle_workflow_error()`, `setup_logging()`, `track_progress()`, `validate_parameters()`, `validate_properties_with_synonyms()`, `validate_data_quality()`, `is_unsuitable()`, `detect_outliers()`, `get_default_configuration()`, `merge_configurations()`.

**External packages**: base `stats` (`qnorm`, `pnorm`, `runif`, `chol`, `sd`, `setNames`, etc.), `parallel` (`detectCores`, `makeCluster`/`stopCluster`/`clusterExport`/`clusterEvalQ`/`parLapply` on Windows, `mclapply` elsewhere). `distributions.R`'s `estimate_correlation_matrix_robust()` in turn depends on `Hmisc::rcorr()`.

**Downstream consumers of this file's output**:
- **GP depth modeling & multivariate adjustment** (`gp-modeling.R`, `multivariate-adjustment.R`) - call `generate_monte_carlo_realizations()` directly per cokey (choosing a single flat correlation matrix via `select_simulation_correlation_matrix()` from a genhz-keyed set), then flatten the resulting `[horizon, property, realization]` array to long format (`result$simulation_data`) for depth-trend and NRCS/local GP adjustment.
- **Statistics & Diagnostics** (`validation-diagnostics.R`) - consumes the long-format simulation data (columns per property plus `simulation_number`, `hzdept_r`, `cokey`, etc.) for quality-assessment reporting, depth-binned summaries, and cross-property/cross-depth diagnostics.

## Data Flow In/Out

**In**:
- `soil_data` - a data frame of SSURGO (or SSURGO-shaped) horizon data with `cokey` and, per requested property, `_l`/`_r`/`_h` percentile-triplet columns; optionally `hzdepb_r` (depth filtering), `hzname` (unsuitability/genhz classification), and `genhz` (explicit genetic-horizon grouping for correlation estimation).
- `properties` - the property name vector to simulate (raw SSURGO names; texture triplets are auto-detected as a composition group when all three are requested).
- Optional `correlation_matrix` (user-supplied, named/dimensioned over `properties`) or auto-estimation via `auto_correlation`/`correlation_fallback` config.
- Optional `observed_data` - field/lab observations (raw samples or parametric summaries) to Bayesian-fuse into the SSURGO-derived priors.
- `simulation_config` - flat or nested tuning overrides (distribution families, depth cutoff, constraint toggles, quality thresholds, parallelism, etc.).

**Out**:
- `simulation_data` - the primary output: a `[horizon, property, realization]` numeric array of constrained, correlated Monte Carlo realizations, dimnamed over the caller-facing `properties` (texture columns already restored from ILR space and guaranteed to sum to 100).
- Supporting metadata: `original_data` (filtered input), `parameters` (fitted per-horizon distributions), `correlation_structure`, `distribution_setup`, `validation` (input/property/output), `diagnostics` (sim-vs-original mean/SD per property), `quality_assessment` (overall + component quality scores), and `metadata` (timing, counts, config, success rate).

## Known Limitations

- `generate_component_values()` (used by `sim_component_compositions()`) is a placeholder: it always draws from a `runif(min, max)` distribution regardless of the requested `distribution_type` - the component-composition pipeline is not yet family-aware in the way the main horizon-property pipeline is.

## Usage Example

```r
mc_results <- generate_monte_carlo_realizations(
  soil_data = soil_data,                 # data frame with cokey + _l/_r/_h columns
  properties = c("sandtotal", "claytotal", "silttotal", "dbovendry"),
  n_realizations = 5000,
  simulation_config = list(               # flat shape is auto-normalized
    max_depth = 200,
    auto_correlation = TRUE,
    correlation_fallback = "kssl_global"
  ),
  parallel = TRUE,
  seed = 12345,
  observed_data = list(
    dbovendry = list(mean = 1.35, sd = 0.08)   # lab-measured bulk density prior update
  )
)

# [horizon, property, realization] array of correlated realizations
str(mc_results$simulation_data)
mc_results$metadata$success_rate
mc_results$quality_assessment$overall_quality_score
```
