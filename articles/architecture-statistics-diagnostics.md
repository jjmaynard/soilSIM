# Statistics & Diagnostics

## Overview

This functional group covers two files — `R/statistics.R` and
`R/validation-diagnostics.R` — that together provide the package’s
analytical and quality-assurance layer. `statistics.R` performs
statistical characterization of processed SSURGO soil property data:
correlation analysis (including compositional/ILR-aware texture
correlations), distribution fitting via `fitdistrplus`, outlier
detection (univariate IQR/Z-score and multivariate Mahalanobis), and
descriptive property statistics with confidence intervals.
`validation-diagnostics.R` performs end-to-end quality assurance across
the *entire* simulation workflow, not just this group’s own output: it
checks Monte Carlo convergence and coverage, correlation-structure
preservation through the simulation pipeline, Gaussian-process (GP)
depth-trend model performance, and soil-science plausibility of final
simulated values, then rolls all of that into a single weighted quality
score, grade, and optional diagnostic report/plots. In the overall
soilSIM pipeline this group sits downstream of data acquisition and
processing: it consumes cleaned/infilled SSURGO data for its own
statistical analysis, and separately consumes the *outputs* of the Monte
Carlo, GP-modeling, and multivariate-adjustment (correlation) groups to
validate the workflow as a whole.

## Core Functions

### Statistical Analysis (`statistics.R`)

#### 1. `analyze_soil_statistics()` — Master Statistical Analysis Entry Point

**Purpose**: Top-level driver that runs the full statistics pipeline
(validation → correlation → distribution fitting → outlier detection →
property statistics → results validation → quality report) over a
processed SSURGO data frame.

**Signature**:

``` r

analyze_soil_statistics(
  processed_data,
  analysis_config = list(),
  correlation_methods = c("pearson", "spearman"),
  distribution_fitting = TRUE,
  outlier_detection = TRUE,
  validate_results = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

**Parameters**: - `processed_data` — processed SSURGO data frame (from
the data-processing group), expected to carry `cokey`, `hzdept_r`, and
numeric `_r` property columns. - `analysis_config` — list of overrides
merged (via
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md))
on top of
[`get_statistical_analysis_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_statistical_analysis_defaults.md). -
`correlation_methods` — character vector of correlation methods to
compute (`"pearson"`, `"spearman"`, and/or `"kendall"`). -
`distribution_fitting` — whether to run Step 5 (distribution
analysis). - `outlier_detection` — whether to run Step 6 (outlier
analysis). - `validate_results` — whether to run Step 8
(statistical-results validation). - `verbose` — enables `DEBUG`-level
logging and initializes logging if not already configured.

**Returns**: A list with `correlation_matrices`,
`distribution_analysis`, `outlier_analysis`, `property_statistics`,
`validation_results`, `analysis_metadata` (timing, config, counts),
`quality_report`, `data_validation`, and `processed_data` (only
populated if `config$return_processed_data` is `TRUE`).

**Behavior**: Runs a fixed 10-step pipeline: (1) validates input data
quality via
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
warning rather than stopping if quality is below `min_quality_score`;
(2) identifies and validates numeric soil properties via
[`identify_numeric_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/identify_numeric_soil_properties.md) +
[`filter_valid_numeric_properties()`](https://jjmaynard.github.io/soilSIM/reference/filter_valid_numeric_properties.md)
(drops all-NA/non-finite columns) +
[`validate_properties_with_synonyms()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties_with_synonyms.md);
(3) handles missing values via
[`handle_missing_values()`](https://jjmaynard.github.io/soilSIM/reference/handle_missing_values.md),
optionally grouped by `genhz`; (4) runs correlation analysis via the
internal
[`run_comprehensive_correlation_analysis_safe()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis_safe.md);
(5) optionally fits distributions via
[`analyze_property_distributions_safe()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions_safe.md)
(a lightweight version that only computes descriptive summary statistics
per property, not actual distribution fits); (6) optionally detects
outliers via
[`detect_comprehensive_outliers_safe()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers_safe.md)
(IQR method only); (7) computes full property statistics via
[`compute_property_statistics()`](https://jjmaynard.github.io/soilSIM/reference/compute_property_statistics.md);
(8) optionally validates the accumulated results via
[`validate_statistical_results_safe()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results_safe.md);
(9) builds a quality report via
[`generate_statistical_quality_report_safe()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report_safe.md);
(10) assembles timing/count metadata. Every analytical step is
individually `tryCatch`-wrapped with
`handle_workflow_error(..., "warn")`, so a failure in any one step
degrades gracefully rather than aborting the whole analysis.

#### 2. `run_comprehensive_correlation_analysis()` — Enhanced Correlation Analysis

**Signature**:
`run_comprehensive_correlation_analysis(data, methods, config, available_properties)`

**Parameters**: `data` (input data frame), `methods` (correlation
methods), `config` (analysis configuration — reads
`config$handle_constant_variables`, `config$stratify_by_horizon`,
`config$include_texture_analysis`), `available_properties` (properties
to correlate).

**Returns**: List with `matrices` (per-method matrix + eigenvalues +
condition number via `calculate_condition_number()`), `stratified` (if
`config$stratify_by_horizon` and a `genhz` column exists, via
[`compute_stratified_correlations()`](https://jjmaynard.github.io/soilSIM/reference/compute_stratified_correlations.md)),
`texture_analysis` (if `config$include_texture_analysis` and ≥2 texture
properties are available, via
[`analyze_texture_correlations()`](https://jjmaynard.github.io/soilSIM/reference/analyze_texture_correlations.md)),
`validation` (via
[`validate_correlation_matrices()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrices.md)),
and `summary`.

**Behavior**: This is a richer, more feature-complete correlation
routine than the internal `_safe` version used by
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)
(see Internal Connections below) — it additionally computes
eigen-decomposition diagnostics per matrix and adds horizon-stratified
and texture-specific correlation blocks. It does not filter to
complete-case rows the way the `_safe` version does, so it assumes
`available_properties` are already appropriately cleaned.

#### 3. `compute_stratified_correlations()`

**Signature**:
`compute_stratified_correlations(data, properties, stratify_by, methods, config)`

**Purpose/Behavior**: Splits `data` by unique non-`NA` values of
`stratify_by` (e.g. `genhz`), skips any group with fewer than
`config$minimum_observations` rows, and computes a
[`safe_correlation()`](https://jjmaynard.github.io/soilSIM/reference/safe_correlation.md)
matrix per method within each remaining group. **Returns** a named list
keyed by group value, each containing one entry per method (`matrix`,
`n_observations`, `group`).

#### 4. `analyze_property_distributions()` — Enhanced Distribution Analysis

**Signature**:
`analyze_property_distributions(data, properties, config)`

**Purpose/Behavior**: For each property with at least
`config$minimum_observations` finite values, calls
[`fit_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/fit_property_distributions.md)
to fit candidate distributions and
[`perform_distribution_tests()`](https://jjmaynard.github.io/soilSIM/reference/perform_distribution_tests.md)
for goodness-of-fit statistics, then summarizes via
`assess_overall_distributions()`. **Returns**
`list(fitted_distributions=, distribution_tests=, summary=)`. Unlike the
internal
[`analyze_property_distributions_safe()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions_safe.md)
used by the master pipeline (which only computes mean/sd/quantiles),
this is the version that performs actual parametric distribution
fitting.

#### 5. `fit_property_distributions()`

**Signature**:
`fit_property_distributions(values, property_name, config)`

**Purpose/Behavior**: Determines candidate distribution families via
[`get_appropriate_distributions()`](https://jjmaynard.github.io/soilSIM/reference/get_appropriate_distributions.md),
fits each with
[`fit_single_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fit_single_distribution.md)
(skipping failures), and ranks the successful fits best-first by AIC via
[`rank_distribution_fits()`](https://jjmaynard.github.io/soilSIM/reference/rank_distribution_fits.md).
**Returns** a named list of fit results (one per successfully-fit
distribution family).

#### 6. `get_appropriate_distributions()`

**Signature**:
`get_appropriate_distributions(property_name, values, config = NULL)`

**Purpose**: Chooses candidate distribution families for a property
using a type-based heuristic — `c("beta", "normal")` for texture
fractions (`sandtotal_r`/`silttotal_r`/`claytotal_r`),
`c("normal", "gamma")` for `ph1to1h2o_r`, and
`c("normal", "lognormal", "gamma", "weibull")` otherwise. If
`config$distribution_methods` is supplied, the result is narrowed to the
intersection of the heuristic and the requested list; if that
intersection is empty, the requested list is used unfiltered (with a
logged warning) rather than silently discarded. `values` is accepted for
interface symmetry with callers but is not currently used by the
heuristic itself. **Returns** a character vector of candidate
distribution names.

#### 7. `detect_comprehensive_outliers()` — Enhanced Outlier Detection

**Signature**: `detect_comprehensive_outliers(data, properties, config)`

**Purpose/Behavior**: For each property with sufficient observations,
runs every method in `config$outlier_methods` (e.g. `"iqr"`, `"zscore"`,
`"modified_zscore"`) via
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md),
using the matching threshold from `config$outlier_thresholds`. If
`length(properties) >= 2` and `config$detect_multivariate_outliers` is
`TRUE`, also runs
[`detect_multivariate_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_multivariate_outliers.md)
(real Mahalanobis-distance detection, guarding a non-positive-definite
covariance matrix with
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html); flags
the outer `1 - config$mahalanobis_alpha` fraction, default 2.5%, via a
chi-squared cutoff). **Returns**
`list(property_outliers=, multivariate_outliers=, summary=)` (summary
via
[`generate_outlier_summary()`](https://jjmaynard.github.io/soilSIM/reference/generate_outlier_summary.md)).

#### 8. `validate_statistical_results()` — Enhanced Results Validation

**Signature**:
`validate_statistical_results(correlation_analysis, distribution_analysis, outlier_analysis, property_statistics, config)`

**Purpose/Behavior**: Validates correlation matrices via
[`validate_correlation_matrices()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrices.md)
(real check, delegates to `distributions.R`’s
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)),
plus pass-through stubs `validate_distribution_analysis()` and
`validate_outlier_analysis()` (both currently always return an empty
warnings list with no real checks). Aggregates into an overall validity
flag and a `validation_score` via `calculate_overall_validation_score()`
(1.0 minus per-error/per-warning penalties). **Returns**
`list(overall_valid=, errors=, warnings=, validation_score=, component_validations=)`.

#### 9. `generate_statistical_quality_report()` — Enhanced Quality Reporting

**Signature**:
`generate_statistical_quality_report(original_data, processed_data, correlation_analysis, distribution_analysis, outlier_analysis, property_statistics, validation_results, data_validation, config)`

**Purpose/Behavior**: Assembles a composite report:
`overall_quality_score` via
[`calculate_overall_statistics_quality_score()`](https://jjmaynard.github.io/soilSIM/reference/calculate_overall_statistics_quality_score.md)
(mean of whichever of the data-quality score, validation score, and
per-component quality assessments are available), `data_quality`
(including
[`calculate_missing_data_improvement()`](https://jjmaynard.github.io/soilSIM/reference/calculate_missing_data_improvement.md)),
`analysis_quality` (via
[`assess_correlation_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_quality.md)
/
[`assess_distribution_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_distribution_quality.md)
/
[`assess_outlier_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_outlier_quality.md)),
`validation_summary`, and `recommendations` (via
`generate_analysis_recommendations()`, currently two hardcoded checks:
data-quality-below-0.8 and presence of correlation warnings).

#### 10. `get_statistical_analysis_defaults()`

**Signature**:
[`get_statistical_analysis_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_statistical_analysis_defaults.md)
(no arguments)

**Purpose/Behavior**: Builds the default statistical-analysis
configuration by taking `get_default_configuration("full")` and merging
in a `statistical_analysis` block covering stratification flags,
missing-value strategy, correlation method/threshold defaults,
distribution-fitting families, outlier method/threshold defaults, and
quality-control thresholds. **Returns** the merged configuration list
(via
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md)).

#### 11. `validate_statistical_config()`

**Signature**: `validate_statistical_config(config)`

**Purpose/Behavior**: Declares a `param_specs` table (types, choices,
min/max lengths, numeric ranges) for the statistical-analysis-specific
config keys (`minimum_observations`, `missing_data_threshold`,
`min_quality_score`, `correlation_methods`, `missing_value_strategy`,
`error_recovery_action`, `outlier_methods`, `distribution_methods`) and
validates via `validate_parameters(..., strict_mode = FALSE)`.
Transparently unwraps a nested `config$statistical_analysis` block if
present. **Returns** the validation result from
[`validate_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_parameters.md).

#### 12. `identify_numeric_soil_properties()`

**Signature**: `identify_numeric_soil_properties(data)`

**Purpose/Behavior**: Finds numeric columns whose name matches `"_r$"`
(representative-value suffix convention), then reorders the result so
that a fixed list of 16 “important” soil properties (`sandtotal_r`,
`silttotal_r`, `claytotal_r`, `dbovendry_r`, `dbthirdbar_r`,
`ph1to1h2o_r`, `cec7_r`, `om_r`, `rfv_r`, `wthirdbar_r`,
`wfifteenbar_r`, `awc_r`, `ksat_r`, `ec_r`, `sar_r`, `esp_r`) come first
(if present), followed by any other matching properties. **Returns** a
character vector of column names.

#### 13. `compute_property_statistics()`

**Signature**: `compute_property_statistics(data, properties, config)`

**Purpose/Behavior**: For each property, computes
n/missing-rate/mean/sd/min/max/range/CV/median/quartiles/IQR on finite
values, adds a 95% t-based confidence interval via
[`calculate_confidence_intervals()`](https://jjmaynard.github.io/soilSIM/reference/calculate_confidence_intervals.md),
and adds skewness/kurtosis
(`calculate_skewness()`/`calculate_kurtosis()`, both hand-rolled
moment-based implementations) plus a Shapiro-Wilk normality p-value
(only computed when n ≤ 5000, since
[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html) is not
defined for larger samples). Properties with zero finite values are
silently skipped. **Returns** a named list keyed by property.

#### 14. `analyze_texture_correlations()`

**Signature**:
`analyze_texture_correlations(data, texture_properties, methods, config)`

**Purpose/Behavior**: Computes raw Pearson/Spearman correlations among
available texture properties (`raw_correlations`), which are known to be
spuriously negative for compositional (simplex-constrained, sum-to-100)
data. When all three of `claytotal_r`/`sandtotal_r`/`silttotal_r` are
present with \>4 complete rows, additionally computes the same
correlations in isometric-log-ratio (ILR) space via `distributions.R`’s
dependency-free
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md),
exposed as `ilr_correlations` — the statistically defensible view for
compositional data. **Returns**
`list(raw_correlations=, ilr_correlations=, n_observations=, note=)`.

**Minor exported helpers in `statistics.R`** (each documented above only
in brief, or not separately elaborated):
[`filter_valid_numeric_properties()`](https://jjmaynard.github.io/soilSIM/reference/filter_valid_numeric_properties.md)
(drops properties with no finite values, non-exported), and internal
`_safe` counterparts
[`run_comprehensive_correlation_analysis_safe()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis_safe.md),
[`analyze_property_distributions_safe()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions_safe.md),
[`detect_comprehensive_outliers_safe()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers_safe.md),
[`validate_statistical_results_safe()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results_safe.md),
[`generate_statistical_quality_report_safe()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report_safe.md)
(all non-exported; these are the lighter-weight versions actually wired
into
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)).
Small numeric utilities `calculate_condition_number()`,
`calculate_skewness()`, `calculate_kurtosis()` are also non-exported.

### Workflow Validation & Diagnostics (`validation-diagnostics.R`)

#### 1. `validate_complete_workflow()` — Master Validation Entry Point

**Purpose**: Top-level driver that assesses an entire soil-simulation
workflow’s output (Monte Carlo, correlation, GP models, soil-science
realism) and rolls it up into one overall quality assessment.

**Signature**:

``` r

validate_complete_workflow(
  workflow_results,
  original_data = NULL,
  validation_config = NULL,
  generate_plots = TRUE,
  output_dir = NULL
)
```

**Parameters**: - `workflow_results` — complete workflow output,
expected to optionally contain `simulation_data`, `integrated_data`,
`gp_models`, `correlation_matrices`, `training_data` keys. -
`original_data` — original SSURGO/NRCS data, used as a comparison
baseline for coverage and depth-trend checks. - `validation_config` —
validation configuration; defaults to
`get_default_configuration("validation")`. - `generate_plots` — whether
to generate diagnostic plots (Step 7 of the pipeline). - `output_dir` —
directory to save diagnostic plot files into (plots are still returned
as objects/paths even if `NULL`, just not written to disk).

**Returns**: A list with `workflow_summary`, `monte_carlo_validation`,
`correlation_validation`, `gp_model_validation`,
`soil_science_validation`, `performance_metrics`, `overall_assessment`
(score/grade/status), `diagnostic_plots`, `recommendations`, and
`validation_metadata` (timestamp, config, R/package versions, duration).

**Behavior**: (1) initializes the result structure via
`initialize_validation_structure()`; (2) extracts and shallow-validates
whichever of
`simulation_data`/`gp_models`/`correlation_matrices`/`training_data` are
present via `extract_and_validate_components()`; (3) runs
`execute_validation_pipeline()`, which in turn conditionally calls
[`validate_monte_carlo_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_monte_carlo_quality.md),
[`validate_correlation_structures()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_structures.md),
[`validate_gp_model_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_model_workflow.md),
and
[`validate_soil_science_realism()`](https://jjmaynard.github.io/soilSIM/reference/validate_soil_science_realism.md)
— each skipped (not failed) if its required component wasn’t found —
then computes
[`calculate_workflow_performance()`](https://jjmaynard.github.io/soilSIM/reference/calculate_workflow_performance.md),
[`assess_workflow_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_workflow_quality.md),
optionally
[`generate_comprehensive_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_comprehensive_diagnostics.md),
and `generate_workflow_recommendations()`; (4) finalizes timing via
`finalize_validation_results()`. Every pipeline step is individually
`tryCatch`-wrapped, so a missing or malformed component degrades that
step to a `validation_skipped`/`validation_failed` entry rather than
aborting the whole call.

#### 2. `generate_validation_report()`

**Signature**:

``` r

generate_validation_report(
  validation_results,
  output_format = "html",
  output_file = NULL,
  include_plots = TRUE
)
```

**Parameters**: `validation_results` (output of
[`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md)),
`output_format` (`"html"`, `"pdf"`, or `"markdown"`), `output_file`
(defaults to a timestamped filename in the current directory),
`include_plots` (whether to embed
`validation_results$diagnostic_plots`).

**Returns**: `TRUE`/`FALSE` success status.

**Behavior**: Validates `output_format`/`include_plots` via
`validate_parameters(strict_mode = TRUE)`, builds report content via
`create_report_content()` (executive summary + detailed results +
plots + recommendations + metadata), then dispatches to
[`generate_html_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_html_report.md),
[`generate_pdf_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_pdf_report.md),
or
[`generate_markdown_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_markdown_report.md).
All three renderers share a dependency-free markdown-lines fallback
([`render_report_as_markdown_lines()`](https://jjmaynard.github.io/soilSIM/reference/render_report_as_markdown_lines.md)
/
[`render_value_as_lines()`](https://jjmaynard.github.io/soilSIM/reference/render_value_as_lines.md))
and only use `rmarkdown`/`pandoc`/`tinytex` (all Suggests-only) when
actually available — HTML falls back to a `<pre>`-wrapped escaped-text
page, PDF falls back to
[`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html) with
paginated monospace text.

#### 3. `assess_workflow_quality()`

**Signature**:
`assess_workflow_quality(validation_results, validation_config)`

**Purpose/Behavior**: Computes a per-component quality score (Monte
Carlo, correlation, GP models, soil science) via
`calculate_component_quality_scores()`, combines them into a single
weighted score via `calculate_weighted_quality_score()` (weights from
`validation_config$quality_weights` or `get_default_quality_weights()`,
an equal 0.25/0.25/0.25/0.25 split by default), assigns a letter-style
grade via `determine_quality_grade()` (`"Excellent"` ≥0.9 down to
`"Poor"` \<0.6), and flags `identify_critical_issues()` (validation
failures, sub-0.6 overall score, data-quality failures). **Returns** a
`create_quality_assessment()` list: `quality_score`, `quality_grade`,
`component_scores`, `critical_issues`, `workflow_status` (`"PASSED"` if
`quality_score >= validation_config$minimum_acceptable_score` — default
0.7 — else `"FAILED"`), `confidence_level` (via
`calculate_confidence_level()`, which rewards a high mean score with low
cross-component variance), and `assessment_metadata`.

#### 4. `validate_monte_carlo_quality()`

**Signature**:
`validate_monte_carlo_quality(monte_carlo_results, original_data = NULL, config = NULL)`

**Purpose/Behavior**: Validates `monte_carlo_results$simulation_data`
via
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)
(requires `cokey`, `simulation_number`), then runs four assessments in
sequence with progress tracking
([`track_progress()`](https://jjmaynard.github.io/soilSIM/reference/track_progress.md)):
[`assess_simulation_convergence()`](https://jjmaynard.github.io/soilSIM/reference/assess_simulation_convergence.md)
(splits `simulation_number` into `config$convergence_criteria$n_batches`
sequential batches — default 5 — and checks whether each numeric
property’s batch mean has stabilized within `tolerance`, default 0.05
relative change),
[`assess_simulation_coverage()`](https://jjmaynard.github.io/soilSIM/reference/assess_simulation_coverage.md),
[`validate_distribution_fidelity()`](https://jjmaynard.github.io/soilSIM/reference/validate_distribution_fidelity.md),
and a `"diagnostics"` step that calls
[`generate_simulation_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_simulation_diagnostics.md).
**Fixed** (previously a bug): that last call used to pass only 2
positional arguments while
[`generate_simulation_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_simulation_diagnostics.md)
(defined in `monte-carlo.R`) takes 5
(`original_data, simulation_results, properties, correlation_config, config`),
so it always raised an argument-mismatch error, silently caught by the
surrounding `tryCatch` into a
`list(validation_failed = TRUE, error = ...)` placeholder instead of
real diagnostics. The call now passes all 5 args by name, sourced from
this function’s own `original_data` parameter, the already-fetched
`simulation_data` local (the constrained
`[horizon, property, realization]` array), and
`monte_carlo_results$metadata$properties`/`monte_carlo_results$correlation_structure`.
**Returns**
`list(convergence_assessment=, coverage_assessment=, distribution_fidelity=, simulation_diagnostics=, data_quality=)`.

#### 5. `assess_simulation_coverage()`

**Signature**:
`assess_simulation_coverage(simulation_data, original_data, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(min_coverage_percentile = 0.95, max_extrapolation_factor = 1.2, min_samples_per_group = 10)`.
Identifies numeric properties (excluding
`simulation_number`/`hzdept_r`/`hzdepb_r`) via
[`validate_properties()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties.md),
then for each: if `original_data` is supplied, computes
[`assess_property_coverage()`](https://jjmaynard.github.io/soilSIM/reference/assess_property_coverage.md)
(overlap between simulated and original ranges, plus
extrapolation-beyond-range detection); always computes
[`assess_distributional_coverage()`](https://jjmaynard.github.io/soilSIM/reference/assess_distributional_coverage.md)
(IQR-based outlier percentage and 5th–95th percentile spread relative to
full range). If `cokey` is present, also computes
[`assess_group_representation()`](https://jjmaynard.github.io/soilSIM/reference/assess_group_representation.md)
(fraction of `cokey` groups meeting `min_samples_per_group`).
**Returns**
`list(parameter_space_coverage=, distributional_coverage=, group_representation=)`.

#### 6. `validate_distribution_fidelity()`

**Signature**:
`validate_distribution_fidelity(simulation_data, simulation_metadata, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(ks_test_alpha = 0.05, moment_tolerance = 0.1, quantile_tolerance = 0.05)`.
For each property named in
`simulation_metadata$distribution_parameters`, draws a reference sample
from the fitted distribution via `distributions.R`’s
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md)
and compares it against the simulated values three ways:
[`test_distribution_fidelity()`](https://jjmaynard.github.io/soilSIM/reference/test_distribution_fidelity.md)
(two-sample Kolmogorov-Smirnov test),
[`compare_distribution_moments()`](https://jjmaynard.github.io/soilSIM/reference/compare_distribution_moments.md)
(relative mean/variance difference),
[`compare_distribution_quantiles()`](https://jjmaynard.github.io/soilSIM/reference/compare_distribution_quantiles.md)
(relative difference at the 10th/50th/90th percentiles). Also runs
IQR-based
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md)
per property. **Returns**
`list(distribution_tests=, moment_comparisons=, quantile_comparisons=, outlier_assessment=)`.

#### 7. `validate_correlation_structures()`

**Signature**:
`validate_correlation_structures(correlation_matrices, simulation_data, config = NULL)`

**Purpose/Behavior**: Validates `simulation_data` via
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
then runs four steps with progress tracking:
[`validate_correlation_matrix_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix_quality.md)
(positive-definiteness + condition number per matrix, via
[`assess_correlation_matrix_properties()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_matrix_properties.md)),
[`validate_correlation_preservation_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_preservation_diagnostics.md),
[`validate_within_depth_correlations()`](https://jjmaynard.github.io/soilSIM/reference/validate_within_depth_correlations.md),
[`assess_cholesky_decomposition()`](https://jjmaynard.github.io/soilSIM/reference/assess_cholesky_decomposition.md).
**Returns**
`list(matrix_quality=, preservation_assessment=, depth_specific_validation=, cholesky_validation=, data_quality=)`.

#### 8. `validate_correlation_preservation_diagnostics()`

**Signature**:
`validate_correlation_preservation_diagnostics(original_correlations, simulation_data, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(max_correlation_difference = 0.1, correlation_rmse_threshold = 0.05, min_samples_for_validation = 30)`.
Computes the overall simulated correlation matrix via
[`safe_correlation()`](https://jjmaynard.github.io/soilSIM/reference/safe_correlation.md)
and compares it against
`original_correlations$global_correlation_matrix` on shared properties
via
[`assess_correlation_differences()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_differences.md)
(real element-wise max/mean/RMSE difference). Also recomputes
correlations within up to 5 distinct `hzdept_r` depth values (each
requiring ≥ `min_samples_for_validation` rows) and compares each via
[`assess_depth_correlation_preservation()`](https://jjmaynard.github.io/soilSIM/reference/assess_depth_correlation_preservation.md).
**Returns**
`list(overall_preservation=, depth_specific_preservation=, property_specific_preservation=)`
(the last is currently never populated).

#### 9. `validate_within_depth_correlations()`

**Signature**:
`validate_within_depth_correlations(simulation_data, criteria = NULL)`

**Purpose/Behavior**: Default criteria includes
`depth_bins = c(0, 15, 30, 60, 100, 200)` and
`min_samples_per_bin = 20`. Bins `hzdept_r` via
[`cut()`](https://rdrr.io/r/base/cut.html), computes a correlation
matrix per bin (via
[`safe_correlation()`](https://jjmaynard.github.io/soilSIM/reference/safe_correlation.md))
that meets the minimum sample size, and records mean absolute
correlation plus
[`assess_correlation_matrix_properties()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_matrix_properties.md)
per bin. Also computes
[`assess_correlation_stability_across_depths()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_stability_across_depths.md)
(how much each pairwise correlation’s value varies across bins).
**Returns**
`list(depth_bin_correlations=, correlation_stability=, depth_trend_correlations=)`
(the last is currently never populated).

#### 10. `assess_cholesky_decomposition()`

**Signature**:
`assess_cholesky_decomposition(correlation_matrices, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(reconstruction_tolerance = 1e-10, condition_number_threshold = 1e12, eigenvalue_threshold = 1e-8)`.
For every matrix entry in `correlation_matrices`, runs
[`assess_single_cholesky_decomposition()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_cholesky_decomposition.md)
— an actual [`chol()`](https://rdrr.io/r/base/chol.html) attempt with
Frobenius-norm reconstruction-error measurement (`||L'L - matrix||`),
graded `"excellent"`/`"good"`/`"poor"` by tolerance. **Returns**
`list(decomposition_quality=, numerical_stability=, reconstruction_accuracy=)`
(the latter two sub-lists are not populated beyond
`decomposition_quality`).

#### 11. `validate_gp_model_workflow()`

**Signature**:
`validate_gp_model_workflow(gp_models, training_data, config = NULL)`

**Purpose/Behavior**: Validates `training_data` (if supplied) via
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
then runs four GP-focused steps with progress tracking:
[`validate_gp_model_performance()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_model_performance.md),
[`validate_gp_predictions()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_predictions.md),
[`assess_depth_trend_realism()`](https://jjmaynard.github.io/soilSIM/reference/assess_depth_trend_realism.md),
[`perform_gp_cross_validation()`](https://jjmaynard.github.io/soilSIM/reference/perform_gp_cross_validation.md).
**Returns**
`list(model_performance=, prediction_quality=, depth_trend_realism=, cross_validation=, training_data_quality=)`.

#### 12. `validate_gp_model_performance()`

**Signature**:
`validate_gp_model_performance(gp_models, training_data, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(max_training_rmse = 0.5, min_r_squared = 0.6, max_condition_number = 1e12)`.
For every property/group combination in `gp_models` whose type is
`"stratified_grouped"`, calls
[`assess_single_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_gp_performance.md),
which delegates training RMSE to `gp-modeling.R`’s
`calculate_model_diagnostics()` and derives R² from that RMSE against
the training data’s own variance. Aggregates via
[`calculate_overall_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/calculate_overall_gp_performance.md)
(mean R²/RMSE across all models). **Returns**
`list(individual_model_performance=, overall_performance=)`.

#### 13. `assess_depth_trend_realism()`

**Signature**: `assess_depth_trend_realism(gp_models, criteria = NULL)`

**Purpose/Behavior**: Default criteria includes
`test_depths = seq(0, 200, by = 5)`,
`realistic_ranges = get_realistic_property_ranges()`,
`monotonicity_tolerance = 0.3`,
`gradient_limits = get_realistic_gradients()`. For each
stratified-grouped GP model, predicts a depth trend via
`gp-modeling.R`’s
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
then judges realism via
[`assess_trend_realism()`](https://jjmaynard.github.io/soilSIM/reference/assess_trend_realism.md)
— which reuses `gp-modeling.R`’s `assess_trend_monotonicity()` and
`assess_realistic_values()`, and additionally counts values outside
`criteria$realistic_ranges[[prop]]` (or non-finite predictions if no
range is configured for that property). Violations are aggregated via
[`summarize_constraint_violations()`](https://jjmaynard.github.io/soilSIM/reference/summarize_constraint_violations.md).
**Returns**
`list(trend_predictions=, realism_assessment=, constraint_violations=)`
(`realism_assessment` is currently never separately populated beyond
what’s nested in `trend_predictions`/`constraint_violations`).

#### 14. `validate_gp_predictions()`

**Signature**: `validate_gp_predictions(gp_models, criteria = NULL)`

**Purpose/Behavior**: Default criteria
`list(test_depths = seq(0, 150, by = 10), uncertainty_threshold = 0.5, smoothness_criteria = 0.1)`.
For each stratified-grouped GP model, calls
[`validate_single_gp_predictions()`](https://jjmaynard.github.io/soilSIM/reference/validate_single_gp_predictions.md),
which predicts a depth trend and measures smoothness as one minus the
normalized mean absolute second difference of the predicted curve (a
discrete-roughness measure). **Returns**
`list(prediction_smoothness=, uncertainty_assessment=, extrapolation_behavior=)`
(only `prediction_smoothness` is currently populated).

#### 15. `validate_soil_science_realism()`

**Signature**:
`validate_soil_science_realism(simulation_data, original_data = NULL, config = NULL)`

**Purpose/Behavior**: Validates `simulation_data` via
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
then runs four steps with progress tracking:
[`assess_property_constraints()`](https://jjmaynard.github.io/soilSIM/reference/assess_property_constraints.md),
[`validate_horizon_characteristics()`](https://jjmaynard.github.io/soilSIM/reference/validate_horizon_characteristics.md)
(surface vs. subsurface range compliance plus depth-transition
smoothness),
[`assess_pedological_relationships()`](https://jjmaynard.github.io/soilSIM/reference/assess_pedological_relationships.md)
(checks clay–CEC and organic-matter–depth correlation direction/strength
when those columns are present),
[`validate_simulation_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/validate_simulation_depth_trends.md)
(correlation between simulated and original per-depth-bin means, when
`original_data` is supplied). **Returns**
`list(property_constraints=, horizon_characteristics=, pedological_relationships=, depth_trend_realism=, data_quality=)`.

#### 16. `assess_property_constraints()`

**Signature**:
`assess_property_constraints(simulation_data, criteria = NULL)`

**Purpose/Behavior**: Defaults to `get_default_property_constraints()`,
which builds `property_ranges` from `get_realistic_property_ranges()`
(hardcoded min/max for clay/sand/silt/pH/bulk density/water
retention/OM/CEC/RFV). For each configured property present in the data,
checks range violations (both directly and via `utils.R`’s
[`validate_numeric_ranges()`](https://jjmaynard.github.io/soilSIM/reference/validate_numeric_ranges.md)),
then calls
[`assess_cross_property_constraints()`](https://jjmaynard.github.io/soilSIM/reference/assess_cross_property_constraints.md)
(texture sum-to-100 check plus any explicit `cross_property_rules`) and
[`detect_distribution_anomalies()`](https://jjmaynard.github.io/soilSIM/reference/detect_distribution_anomalies.md)
(IQR-outlier counts per numeric column, via
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md)).
**Returns**
`list(range_violations=, cross_property_violations=, distribution_anomalies=)`.

**Minor exported helpers grouped by theme** (each is a real, working
assessment function but narrower in scope than the 16 above — none
carries its own `@section Known limitation`): - *Monte Carlo internals*:
[`assess_simulation_convergence()`](https://jjmaynard.github.io/soilSIM/reference/assess_simulation_convergence.md),
[`assess_property_coverage()`](https://jjmaynard.github.io/soilSIM/reference/assess_property_coverage.md),
[`assess_distributional_coverage()`](https://jjmaynard.github.io/soilSIM/reference/assess_distributional_coverage.md),
[`assess_group_representation()`](https://jjmaynard.github.io/soilSIM/reference/assess_group_representation.md),
[`test_distribution_fidelity()`](https://jjmaynard.github.io/soilSIM/reference/test_distribution_fidelity.md),
[`compare_distribution_moments()`](https://jjmaynard.github.io/soilSIM/reference/compare_distribution_moments.md),
[`compare_distribution_quantiles()`](https://jjmaynard.github.io/soilSIM/reference/compare_distribution_quantiles.md). -
*Correlation internals*:
[`validate_correlation_matrix_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix_quality.md),
[`assess_correlation_differences()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_differences.md),
[`assess_depth_correlation_preservation()`](https://jjmaynard.github.io/soilSIM/reference/assess_depth_correlation_preservation.md),
[`assess_correlation_matrix_properties()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_matrix_properties.md),
[`assess_correlation_stability_across_depths()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_stability_across_depths.md),
[`assess_single_cholesky_decomposition()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_cholesky_decomposition.md). -
*GP internals*:
[`assess_single_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_gp_performance.md),
[`calculate_overall_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/calculate_overall_gp_performance.md),
[`assess_trend_realism()`](https://jjmaynard.github.io/soilSIM/reference/assess_trend_realism.md),
[`summarize_constraint_violations()`](https://jjmaynard.github.io/soilSIM/reference/summarize_constraint_violations.md),
[`validate_single_gp_predictions()`](https://jjmaynard.github.io/soilSIM/reference/validate_single_gp_predictions.md),
[`perform_gp_cross_validation()`](https://jjmaynard.github.io/soilSIM/reference/perform_gp_cross_validation.md). -
*Soil-science internals*: `get_default_property_constraints()`,
[`assess_cross_property_constraints()`](https://jjmaynard.github.io/soilSIM/reference/assess_cross_property_constraints.md),
[`detect_distribution_anomalies()`](https://jjmaynard.github.io/soilSIM/reference/detect_distribution_anomalies.md),
[`validate_horizon_characteristics()`](https://jjmaynard.github.io/soilSIM/reference/validate_horizon_characteristics.md),
[`assess_pedological_relationships()`](https://jjmaynard.github.io/soilSIM/reference/assess_pedological_relationships.md),
[`validate_simulation_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/validate_simulation_depth_trends.md). -
*Scoring/reporting internals*:
[`calculate_workflow_performance()`](https://jjmaynard.github.io/soilSIM/reference/calculate_workflow_performance.md),
`calculate_monte_carlo_score()`, `calculate_correlation_score()`,
`calculate_gp_model_score()`, `calculate_soil_science_score()`,
`calculate_component_quality_scores()`,
`calculate_weighted_quality_score()`, `determine_quality_grade()`,
`identify_critical_issues()`, `create_quality_assessment()`,
`calculate_confidence_level()`. - *Report rendering/plotting internals*
(non-exported): `create_report_content()`,
[`render_report_as_markdown_lines()`](https://jjmaynard.github.io/soilSIM/reference/render_report_as_markdown_lines.md),
[`render_value_as_lines()`](https://jjmaynard.github.io/soilSIM/reference/render_value_as_lines.md),
[`escape_html_text()`](https://jjmaynard.github.io/soilSIM/reference/escape_html_text.md),
[`generate_html_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_html_report.md),
[`generate_markdown_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_markdown_report.md),
[`generate_pdf_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_pdf_report.md),
[`theme_soil_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/theme_soil_diagnostics.md),
[`save_diagnostic_plot()`](https://jjmaynard.github.io/soilSIM/reference/save_diagnostic_plot.md),
[`generate_comprehensive_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_comprehensive_diagnostics.md),
`generate_workflow_recommendations()`, `create_executive_summary()`,
`create_detailed_results()`. - *Static reference data* (non-exported):
`get_default_quality_weights()`, `get_realistic_property_ranges()`,
`get_realistic_gradients()`.

## Internal Connections

    statistics.R
    ========================================================================
    analyze_soil_statistics()                          [MASTER — statistics.R]
    ├── validate_data_quality()                         [utils.R]
    ├── identify_numeric_soil_properties()
    ├── filter_valid_numeric_properties()
    ├── validate_properties_with_synonyms()              [utils.R]
    ├── handle_missing_values()                          [utils.R]
    ├── run_comprehensive_correlation_analysis_safe()
    │   └── safe_correlation()                           [utils.R]
    ├── analyze_property_distributions_safe()            (descriptive only —
    │                                                       does NOT call fit_property_distributions())
    ├── detect_comprehensive_outliers_safe()
    │   └── detect_outliers()                            [utils.R]
    ├── compute_property_statistics()
    │   ├── calculate_confidence_intervals()             [utils.R]
    │   ├── calculate_skewness() / calculate_kurtosis()
    │   └── shapiro.test()                               [stats]
    ├── validate_statistical_results_safe()
    └── generate_statistical_quality_report_safe()

      NOTE: the "enhanced" exported functions below form a separate,
      richer call chain that is NOT wired into analyze_soil_statistics() —
      they are available for direct use but the master pipeline uses the
      lighter "_safe" internal versions above instead.

    run_comprehensive_correlation_analysis()             [exported, standalone]
    ├── safe_correlation()                               [utils.R]
    ├── compute_stratified_correlations()
    │   └── safe_correlation()                           [utils.R]
    ├── analyze_texture_correlations()
    │   └── ilr_forward()                                [distributions.R]
    └── validate_correlation_matrices()
        └── validate_correlation_matrix()                [distributions.R]

    analyze_property_distributions()                     [exported, standalone]
    ├── fit_property_distributions()
    │   ├── get_appropriate_distributions()
    │   ├── fit_single_distribution()
    │   │   └── fitdistrplus::fitdist()                  [fitdistrplus]
    │   └── rank_distribution_fits()
    └── perform_distribution_tests()
        ├── fit_single_distribution()  (per candidate)
        └── fitdistrplus::gofstat()                       [fitdistrplus]

    detect_comprehensive_outliers()                      [exported, standalone]
    ├── detect_outliers()                                [utils.R]  (per method)
    ├── detect_multivariate_outliers()
    │   └── Matrix::nearPD() / stats::mahalanobis()       [Matrix, stats]
    └── generate_outlier_summary()

    validate_statistical_results()                       [exported, standalone]
    ├── validate_correlation_matrices()
    ├── validate_distribution_analysis()   (stub)
    ├── validate_outlier_analysis()        (stub)
    └── calculate_overall_validation_score()

    generate_statistical_quality_report()                [exported, standalone]
    ├── calculate_overall_statistics_quality_score()
    ├── calculate_missing_data_improvement()
    ├── assess_correlation_quality()
    ├── assess_distribution_quality()
    ├── assess_outlier_quality()
    └── generate_analysis_recommendations()


    validation-diagnostics.R
    ========================================================================
    validate_complete_workflow()                          [MASTER — validation-diagnostics.R]
    ├── initialize_validation_structure()
    ├── extract_and_validate_components()
    │   └── validate_data_quality()                       [utils.R]
    ├── execute_validation_pipeline()
    │   ├── validate_monte_carlo_quality()          (if simulation_data found)
    │   │   ├── assess_simulation_convergence()
    │   │   ├── assess_simulation_coverage()
    │   │   │   ├── assess_property_coverage()
    │   │   │   ├── assess_distributional_coverage()
    │   │   │   │   └── detect_outliers()                 [utils.R]
    │   │   │   └── assess_group_representation()
    │   │   ├── validate_distribution_fidelity()
    │   │   │   ├── test_distribution_fidelity()
    │   │   │   ├── compare_distribution_moments()
    │   │   │   ├── compare_distribution_quantiles()
    │   │   │   │    (all three via quantile_from_fit() [distributions.R])
    │   │   │   └── detect_outliers()                     [utils.R]
    │   │   └── generate_simulation_diagnostics()          [monte-carlo.R — arg-count mismatch, always caught & warned]
    │   ├── validate_correlation_structures()       (if correlation_matrices found)
    │   │   ├── validate_correlation_matrix_quality()
    │   │   │   └── assess_correlation_matrix_properties()
    │   │   ├── validate_correlation_preservation_diagnostics()
    │   │   │   ├── safe_correlation()                     [utils.R]
    │   │   │   ├── assess_correlation_differences()
    │   │   │   └── assess_depth_correlation_preservation()
    │   │   ├── validate_within_depth_correlations()
    │   │   │   ├── safe_correlation()                     [utils.R]
    │   │   │   ├── assess_correlation_matrix_properties()
    │   │   │   └── assess_correlation_stability_across_depths()
    │   │   └── assess_cholesky_decomposition()
    │   │       └── assess_single_cholesky_decomposition()
    │   ├── validate_gp_model_workflow()            (if gp_models found)
    │   │   ├── validate_gp_model_performance()
    │   │   │   ├── assess_single_gp_performance()
    │   │   │   │   └── calculate_model_diagnostics()      [gp-modeling.R]
    │   │   │   └── calculate_overall_gp_performance()
    │   │   ├── validate_gp_predictions()
    │   │   │   └── validate_single_gp_predictions()
    │   │   │       └── predict_gp_depth_trends()          [gp-modeling.R]
    │   │   ├── assess_depth_trend_realism()
    │   │   │   └── assess_trend_realism()
    │   │   │       ├── assess_trend_monotonicity()        [gp-modeling.R]
    │   │   │       └── assess_realistic_values()          [gp-modeling.R]
    │   │   └── perform_gp_cross_validation()
    │   │       └── k_fold_gp_cv()                         [gp-modeling.R]
    │   ├── validate_soil_science_realism()         (if final_data found)
    │   │   ├── assess_property_constraints()
    │   │   │   ├── get_default_property_constraints()
    │   │   │   ├── assess_cross_property_constraints()
    │   │   │   └── detect_distribution_anomalies()
    │   │   │       └── detect_outliers()                  [utils.R]
    │   │   ├── validate_horizon_characteristics()
    │   │   ├── assess_pedological_relationships()
    │   │   └── validate_simulation_depth_trends()
    │   ├── calculate_workflow_performance()
    │   ├── assess_workflow_quality()
    │   │   ├── calculate_component_quality_scores()
    │   │   │   ├── calculate_monte_carlo_score()
    │   │   │   ├── calculate_correlation_score()
    │   │   │   ├── calculate_gp_model_score()
    │   │   │   └── calculate_soil_science_score()
    │   │   ├── calculate_weighted_quality_score()
    │   │   ├── determine_quality_grade()
    │   │   ├── identify_critical_issues()
    │   │   └── create_quality_assessment()
    │   │       └── calculate_confidence_level()
    │   ├── generate_comprehensive_diagnostics()    (if generate_plots)
    │   │   └── ggplot2 / corrplot (Suggests-only, degrades to empty list)
    │   └── generate_workflow_recommendations()
    └── finalize_validation_results()

    generate_validation_report()                          [standalone]
    ├── validate_parameters()                              [utils.R]
    ├── create_report_content()
    │   ├── create_executive_summary()
    │   └── create_detailed_results()
    └── generate_html_report() / generate_pdf_report() / generate_markdown_report()
        └── render_report_as_markdown_lines()
            └── render_value_as_lines()

## Dependencies

**From within soilSIM**: - `R/utils.R` (“Module 0” in legacy terms) —
the single heaviest dependency for both files:
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md),
[`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md),
[`setup_logging()`](https://jjmaynard.github.io/soilSIM/reference/setup_logging.md),
[`get_default_configuration()`](https://jjmaynard.github.io/soilSIM/reference/get_default_configuration.md),
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md),
[`validate_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_parameters.md),
[`validate_properties()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties.md),
[`validate_properties_with_synonyms()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties_with_synonyms.md),
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
[`validate_numeric_ranges()`](https://jjmaynard.github.io/soilSIM/reference/validate_numeric_ranges.md),
[`handle_missing_values()`](https://jjmaynard.github.io/soilSIM/reference/handle_missing_values.md),
[`safe_correlation()`](https://jjmaynard.github.io/soilSIM/reference/safe_correlation.md),
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md),
[`calculate_confidence_intervals()`](https://jjmaynard.github.io/soilSIM/reference/calculate_confidence_intervals.md),
[`track_progress()`](https://jjmaynard.github.io/soilSIM/reference/track_progress.md),
`get_default_quality_thresholds()`, `get_package_versions()`. -
`R/distributions.R` —
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)
(isometric log-ratio transform for compositional texture correlation),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
(singular; the real per-matrix check backing
[`validate_correlation_matrices()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrices.md),
plural, in `statistics.R`),
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md)
(theoretical quantiles from a fitted distribution, used throughout the
Monte Carlo fidelity checks). - `R/gp-modeling.R` —
`calculate_model_diagnostics()`,
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
`assess_trend_monotonicity()`, `assess_realistic_values()`,
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
— all consumed only by `validation-diagnostics.R`’s GP-validation
functions. - `R/monte-carlo.R` —
[`generate_simulation_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_simulation_diagnostics.md),
called from
[`validate_monte_carlo_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_monte_carlo_quality.md)
(argument-count mismatch fixed — see
[`validate_monte_carlo_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_monte_carlo_quality.md)
above).

**External packages**: - `fitdistrplus` — `fitdist()` and `gofstat()`,
the entire basis of
[`fit_single_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fit_single_distribution.md)
/
[`perform_distribution_tests()`](https://jjmaynard.github.io/soilSIM/reference/perform_distribution_tests.md).
Required ([`stop()`](https://rdrr.io/r/base/stop.html)s with an install
hint if missing) rather than optional, despite being invoked lazily via
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html). - `stats` —
correlation, [`cor()`](https://rdrr.io/r/stats/cor.html),
[`cov()`](https://rdrr.io/r/stats/cor.html),
[`mahalanobis()`](https://rdrr.io/r/stats/mahalanobis.html),
[`chol()`](https://rdrr.io/r/base/chol.html),
[`eigen()`](https://rdrr.io/r/base/eigen.html),
[`quantile()`](https://rspatial.github.io/terra/reference/quantile.html),
[`sd()`](https://rdrr.io/r/stats/sd.html),
[`var()`](https://rdrr.io/r/stats/cor.html),
[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html),
[`ks.test()`](https://rdrr.io/r/stats/ks.test.html),
[`qchisq()`](https://rdrr.io/r/stats/Chisquare.html),
[`weighted.mean()`](https://rspatial.github.io/terra/reference/weighted.mean.html),
[`complete.cases()`](https://rdrr.io/r/stats/complete.cases.html). -
`Matrix` — `nearPD()`, used to repair a non-positive-definite covariance
matrix before Mahalanobis distance computation in
[`detect_multivariate_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_multivariate_outliers.md). -
`Hmisc` (declared package Import, used indirectly via `utils.R`’s
statistical helpers such as
[`calculate_confidence_intervals()`](https://jjmaynard.github.io/soilSIM/reference/calculate_confidence_intervals.md)). -
`dplyr` / `tidyr` —
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
(depth/bin subsetting in `validation-diagnostics.R`),
[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
/
[`dplyr::everything()`](https://tidyselect.r-lib.org/reference/everything.html)
(Monte Carlo diagnostic plot reshaping). - `ggplot2`, `corrplot`
(Suggests-only) — diagnostic plot rendering in
[`generate_comprehensive_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_comprehensive_diagnostics.md);
every plot block checks
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) first and
degrades to an empty list rather than erroring if unavailable. -
`rmarkdown`, `tinytex` (Suggests-only) — richer HTML/PDF report
rendering in
[`generate_html_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_html_report.md)/[`generate_pdf_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_pdf_report.md),
with a dependency-free fallback in both cases. - `grDevices`, `graphics`
— PDF/PNG diagnostic plot fallbacks.

**Upstream/downstream soilSIM functional groups**: - **Upstream (feeds
this group)**: Data Acquisition & Processing (`ssurgo-acquisition.R`,
`ssurgo-processing.R`, `data-infilling.R`) supplies the `processed_data`
consumed by
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md). -
**Downstream/consumers (validated by this group, not fed by it)**:
[`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md)
and its sub-validators are the quality gate for the Monte Carlo
simulation group (`monte-carlo.R`), the GP depth-trend modeling group
(`gp-modeling.R`), and the multivariate/correlation adjustment group
(`multivariate-adjustment.R`) — this group consumes their *output*
(`simulation_data`, `gp_models`, `correlation_matrices`) purely for
assessment; it does not feed data back into them.

## Data Flow In/Out

**In**: -
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)
— a processed SSURGO data frame (post-infilling), with `cokey`,
`hzdept_r`, and numeric `_r` property columns. -
[`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md)
— a `workflow_results` list assembled from one or more upstream stages:
`simulation_data`/`integrated_data` (Monte Carlo or fusion output),
`gp_models` (from GP modeling), `correlation_matrices` (from
correlation/multivariate adjustment), `training_data`; plus optionally
`original_data` for before/after comparison.

**Out**: -
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)
— a nested statistics report: correlation matrices (global + optionally
stratified/texture-specific), distribution summaries, outlier
flags/counts, per-property descriptive statistics with confidence
intervals, a results-validation block, and an overall quality
report/score. -
[`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md)
— a comprehensive validation report object: per-subsystem validation
blocks (Monte Carlo, correlation, GP, soil science), workflow
performance metrics, an overall weighted quality score/grade/pass-fail
status, optional diagnostic plot objects/file paths, and a list of
textual recommendations — suitable for passing directly into
[`generate_validation_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_validation_report.md)
to render as HTML/PDF/Markdown.

## Usage Example

``` r

# Statistical characterization of processed SSURGO data
stats_result <- analyze_soil_statistics(
  processed_data,
  analysis_config = list(
    stratify_by_horizon = TRUE,
    include_texture_analysis = TRUE,
    correlation_methods = c("pearson", "spearman")
  )
)
stats_result$quality_report$overall_quality_score

# End-to-end validation of a completed simulation workflow
validation_result <- validate_complete_workflow(
  workflow_results = list(
    simulation_data = monte_carlo_output$simulation_data,
    gp_models = gp_models,
    correlation_matrices = correlation_matrices
  ),
  original_data = processed_data,
  generate_plots = FALSE
)

validation_result$overall_assessment$quality_grade   # e.g. "Good"
validation_result$overall_assessment$workflow_status  # "PASSED" or "FAILED"

# Render a report from the validation results
generate_validation_report(validation_result, output_format = "markdown",
                            output_file = "validation_report.md")
```
