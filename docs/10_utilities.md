# Utilities

## Overview

`R/utils.R` is the shared-utility foundation for the entire `soilSIM` package. It provides data validation and quality control, logging, configuration management, caching/backup helpers, WKT geometry validation, and property-name/synonym resolution, so that no other functional group in the package (SSURGO acquisition/processing, data infilling, Monte Carlo simulation, GP modeling, multivariate adjustment, statistics, validation/diagnostics) needs to duplicate this logic. It is a **leaf dependency** in the ordinary sense - it sources nothing else in `soilSIM` and is depended on, directly or indirectly, by every other functional group - with one narrow, source-verified exception documented under [Known Limitations](#known-limitations): `get_predefined_properties("ssurgo")` calls `create_ssurgo_property_lookup_working()`, which actually lives in `R/ssurgo-acquisition.R`, not here.

## Core Functions

### Validation & QC

#### `is_unsuitable(data, hzname_col = "hzname", strict_mode = TRUE, custom_exclusions = NULL)`

**Purpose**: Flags soil horizons that should be excluded from property modeling/imputation because they aren't representative soil material (bedrock, cemented pans, organic horizons, restrictive layers, etc.).

**Parameters**:
- `data` - Soil horizon data frame with horizon designation information.
- `hzname_col` - Column name containing horizon designations (default `"hzname"`). If this column is absent, the function logs a WARN and returns all-`FALSE` (nothing flagged) rather than erroring.
- `strict_mode` - Whether to apply additional strict-mode criteria (default `TRUE`, note this differs from the legacy `mod00` default of `FALSE`).
- `custom_exclusions` - Optional additional regex patterns (matched against raw `hzname` values) to flag as unsuitable.

**Returns**: A logical vector, one element per row of `data`, `TRUE` where the horizon is considered unsuitable.

**Algorithm/behavior**: Horizon names are coerced to character and NA-filled to `""`. For each row, `hzname` (and `desgnmaster` if present) are upper-cased/trimmed and checked, in priority order, against: (1) `desgnmaster` in `c("R","O")`; (2) a regex for bedrock/cemented-pan codes (`R`, `Cr`, `Cd`, `Cx`, with an optional leading numeric prefix); (3) restrictive-layer descriptor words (`fragipan`, `duripan`, `petrocalcic`, `petrogypsic`, `cemented`, `indurated`); (4) a trailing cemented-suffix `"m"`; (5) an `O` organic-horizon pattern. If the richer columns `is_potentially_restrictive`, `missing_required_properties`, `hzdept_r`, and `restriction_top` are all present, horizons at or below a restrictive contact (`hzdept_r >= restriction_top`) that are also flagged potentially-restrictive and missing required properties are OR'd into the result. If an `organic_horizon` column is present, it is OR'd in directly. In `strict_mode`, horizons deeper than 250 cm are flagged, and a second pass applies generic transitional-horizon patterns (`^C$`, `^C[0-9]`, `BC`, `CB`, `2C`, `3C`, `4C`) case-sensitively against the raw (non-uppercased) `hznames`. Any `custom_exclusions` regexes are applied last, also case-sensitively against raw `hznames`. Final counts and (if <= 10 distinct) the actual unsuitable horizon-name values are logged via `log_message()`. Note the default `strict_mode = TRUE` here differs from the legacy `mod00_soil_utils.R` prototype's default of `FALSE`, and the strict-mode pattern pass and custom-exclusion pass both match against the **raw, non-uppercased** `hznames` vector rather than the upper-cased `hz` used earlier in the loop - a mixed-case horizon name (e.g. `"bc"`) will not match the strict-mode `"BC"` pattern.

#### `validate_data_quality(data, required_columns = character(0), numeric_columns = character(0), quality_thresholds = NULL)`

**Purpose**: Composite data-quality assessment - missing-data percentage, outlier prevalence (informational only), and consistency checks - rolled into an overall quality score and letter grade.

**Parameters**:
- `data` - Input data frame to assess.
- `required_columns` - Columns whose presence is required; missing ones lower the consistency score.
- `numeric_columns` - Columns expected to be numeric; drives the numeric-completeness and outlier sub-assessments.
- `quality_thresholds` - Optional list overriding `max_missing_percentage` (default 0.2), `min_completeness` (default 0.8), `max_outlier_percentage` (default 0.05), `min_numeric_columns` (default 1); any thresholds omitted from a supplied list are backfilled with these defaults via `%||%`.

**Returns**: A list with `overall_quality` (`score`, `grade`, `validation_passed`, `component_scores`), `column_assessments` (unused/always empty in the current implementation), `missing_data_summary`, `numeric_data_summary`, `outlier_summary`, and `consistency_summary`.

**Algorithm/behavior**: Each of the four assessment steps (missing data, numeric validation, outliers, consistency) is wrapped in its own `tryCatch()` with a safe fallback so one failing step cannot abort the whole assessment. Step 1 computes per-column `mean(is.na(col))` (treating a zero-length column as 100% missing) and flags columns exceeding `max_missing_percentage`. Step 2, only if `numeric_columns` is non-empty, computes per-column total/finite/missing/infinite counts and completeness (`finite/total`) for each existing numeric column. Step 3 runs a simple 1.5x-IQR outlier count *inlined directly in this function* (it does **not** call `detect_outliers()`) across the existing numeric columns, requiring more than 4 finite values per column. Step 4 checks for missing `required_columns` and counts exact duplicate rows via `duplicated(data)`. The overall score is a weighted blend: `missing_score * 0.4 + numeric_score * 0.4 + consistency_score * 0.2`, where `missing_score = 1 - avg_missing_percentage` and `numeric_score` defaults to `0.8` if no numeric summary was computed. The blended score maps to a grade: `>= 0.9` Excellent, `>= 0.8` Good, `>= 0.7` Acceptable, `>= 0.6` Needs Improvement, else Poor; `validation_passed` is `score >= quality_thresholds$min_completeness` (default 0.7... note the default constant used here is `0.7`, not the `min_completeness` default of `0.8` shown above - the inline fallback literal at the pass/fail check differs from the value written into `quality_thresholds$min_completeness` when thresholds are auto-generated).

#### `check_required_columns(data, column_specifications, strict_mode = TRUE)`

**Purpose**: Validates that a set of required columns exist and, optionally, that each has the expected R class.

**Parameters**: `data` - input data frame. `column_specifications` - named list keyed by column name, each entry optionally holding `$type` (an expected `class()[1]` string). `strict_mode` - when `TRUE`, type mismatches are recorded (checked only if `strict_mode` is `TRUE`).

**Returns**: A list with `missing_columns` (character vector), `type_mismatches` (character vector of `"col (expected: X, actual: Y)"` strings), and `validation_passed` (logical).

**Algorithm/behavior**: Iterates the names of `column_specifications`; for each, checks existence in `data` first (recording a miss and skipping further checks for that column via `next` if absent), then compares `class(data[[col_name]])[1]` against `col_spec$type` when both a type is specified and `strict_mode` is on. This function performs no logging and is not currently called elsewhere in the package (see Known Limitations).

#### `validate_numeric_ranges(data, property_ranges, action = "warn")`

**Purpose**: Validates that numeric property columns fall within realistic min/max ranges, with a choice of remediation.

**Parameters**: `data` - input data frame. `property_ranges` - named list keyed by column name, each value a list with `$min`/`$max`. `action` - `"warn"` (log and continue), `"error"` (`stop()` on first violating property), or `"clip"` (clamp out-of-range values in place with `pmax(pmin(...))`).

**Returns**: A list with `range_violations` (named by property, each entry `n_violations`, `violation_indices`, `out_of_range_values`) and `corrected_data` (equal to `data` unless `action = "clip"` was applied).

**Algorithm/behavior**: For each property present in both `data` and `property_ranges`, and numeric, it computes `which(prop_data < min | prop_data > max)`; if any violations exist they are recorded, then handled per `action`. `"error"` stops immediately (later properties are never checked). `"clip"` mutates `corrected_data` in place for that column's violating rows only. `"warn"` just logs.

#### WKT geometry validation

The package ships two independent WKT-validation code paths that do **not** call each other: a monolithic function (`validate_wkt_geometry()`) that is the one actually wired into the rest of `soilSIM`, and a set of granular single-purpose validators that exist but are not invoked from `validate_wkt_geometry()` or from anywhere else in the package (see Known Limitations).

##### `validate_wkt_geometry(wkt_string, crs = "epsg:4326", validation_context = "geographic", area_limits = list(min = 0.0001, max = 100, max_action = "warn"), complexity_limits = list(max_vertices = 1000, max_parts = 100), bounds_check = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90), strict_mode = TRUE, return_geometry = FALSE)`

**Purpose**: The package's actual, load-bearing WKT validator used before an AOI-based SSURGO query is attempted - checks geometry validity, extracts bounding box/area, and screens against configurable area, complexity, and coordinate-bounds limits.

**Parameters**: `wkt_string` - the WKT string to validate. `crs` - CRS string passed to `terra::vect()`. `validation_context` - descriptive label only (not dispatched on internally beyond the defaults callers pass in). `area_limits` - `min`, `max`, `max_action` (`"warn"` downgrades an over-max area to a warning; anything else makes it a hard error). `complexity_limits` - currently accepted as a parameter but **not read anywhere in the function body** (see Known Limitations). `bounds_check` - `xmin`/`xmax`/`ymin`/`ymax` for the bounding-box check. `strict_mode` - accepted but not read in the function body. `return_geometry` - accepted but not read in the function body (the parsed `terra` geometry is never attached to the result even when `TRUE`).

**Returns**: A list with `valid` (logical), `errors`/`warnings` (character vectors), and `geometry_stats` (`bbox`, `area_deg2`, `area_value`, `area_units`, `geometry_type`, `n_vertices` when computable).

**Algorithm/behavior**: Requires `terra` (`requireNamespace()` guard, `stop()`s if absent), then parses `wkt_string` with `terra::vect(wkt_string, crs = crs)` inside a `tryCatch()`. If `terra::is.valid(geom)` is `FALSE`, records an error and returns immediately. Otherwise it extracts coordinates via `terra::crds(geom)`; if that fails or returns a non-matrix, it falls back to a hard-coded default `geometry_stats` (global bbox, area 1.0 deg^2, `"POLYGON"`, 4 vertices) rather than failing. When coordinates are available, it computes the bbox directly from the coordinate matrix's min/max, computes area via `terra::expanse(geom, unit = "m")`, and converts to square-degrees using a fixed approximation of 111,000 m per degree (`area_deg2 <- area_value / (111000 * 111000)`) - this is a rough approximation, not a proper geodesic calculation, and is least accurate away from the equator. It then checks the bbox against `bounds_check` (warning only, never an error) and the area against `area_limits` (below-min is always a warning; above-max is a warning or a hard error depending on `area_limits$max_action`). Any error during the whole `tryCatch` block is caught, logged, and papered over with the same hard-coded default `geometry_stats` used in the coordinate-extraction fallback, so a parse failure still produces a structurally complete (if fictitious) result rather than propagating the exception - callers must inspect `errors`/`warnings`, not rely on an exception, to detect failure.

##### `validate_wkt_string(wkt_string)`

**Purpose**: Cheap, non-`terra`, pure string-level sanity check of a WKT string (no geometry parsing).

**Returns**: List with `valid`, `errors`, `warnings`.

**Algorithm/behavior**: Rejects `NULL`/missing/empty input outright. Rejects non-character input. If a vector of length > 1 is passed, warns and uses only the first element (this warning is generated but the truncation is not actually applied to the string tested further below in the same call, since the local reassignment happens only within that `if` block - the function proceeds to check the *original* `wkt_trimmed` derived from the truncated value, which is consistent). Checks that the trimmed string starts with one of the seven standard geometry-type keywords (`POINT`, `LINESTRING`, `POLYGON`, `MULTIPOINT`, `MULTILINESTRING`, `MULTIPOLYGON`, `GEOMETRYCOLLECTION`) followed by `(` - a miss is a warning, not an error. Checks parenthesis balance by comparing counts of `(` vs `)` via `nchar()` diffing - imbalance is an error. Flags (warning only) WKT strings longer than 100,000 characters. Not called by `validate_wkt_geometry()` or anywhere else in the package.

##### `parse_wkt_geometry(wkt_string, crs)`

**Purpose**: A standalone `terra`-based WKT parser that returns the parsed geometry object plus richer per-geometry statistics (distinct from, and unused by, `validate_wkt_geometry()`).

**Returns**: List with `valid`, `errors`, `warnings`, `geometry` (the `terra` `SpatVector`, always populated on success unlike `validate_wkt_geometry()`), and `stats` (`geometry_type`, `n_geometries`, `n_vertices`, `bbox`, `crs`, `is_valid`, `area_value`, `area_units`).

**Algorithm/behavior**: Parses with `terra::vect(wkt_string, crs = crs)`. Computes vertex count as `sum(terra::nrow(terra::geom(geom)))` and part count as `terra::nrow(geom)`. Area/units branch on `terra::is.lonlat(geom)`: geographic geometries get a crude `(xmax-xmin)*(ymax-ymin)` bbox-area approximation in `deg^2` (not a true polygon area), projected geometries get `terra::expanse(geom, unit = "m")` in `m^2`. On parse failure, inspects the error message text for the substrings `"GDAL"` and `"CRS"` to append more specific guidance to `errors`/`warnings`. Not called by `validate_wkt_geometry()` or elsewhere in the package.

##### `validate_geometry_validity(geom, strict_mode = TRUE)`, `validate_geometry_area(geom, area_limits, context = "general")`, `validate_geometry_complexity(geom, complexity_limits)`, `validate_coordinate_bounds(geom, bounds_check, context = "general")`, `validate_geographic_context(geom, strict_mode = TRUE)`, `validate_projected_context(geom, crs, strict_mode = TRUE)`, `get_validation_defaults(context, crs)`, `calculate_complexity_score(n_vertices, n_parts, geometry_type)`

**Purpose**: A granular, single-responsibility family of geometry checks that together duplicate (with somewhat more nuance) what `validate_wkt_geometry()` does inline. Each takes an already-parsed `terra` geometry (from `parse_wkt_geometry()`, in principle) rather than a raw WKT string.

**Algorithm/behavior (brief, by function)**: `validate_geometry_validity()` calls `terra::is.valid()`; on failure, downgrades to a warning unless `strict_mode`, and additionally attempts `terra::makeValid()` to see if the geometry is auto-fixable, noting that in a warning. `validate_geometry_area()` computes area the same lonlat-branching way as `parse_wkt_geometry()` and compares against `area_limits$min`/`max`, each independently escalatable to an error via `area_limits$min_action`/`max_action == "error"`. `validate_geometry_complexity()` computes vertex/part counts and calls `calculate_complexity_score()`, warning (never erroring) if either exceeds `complexity_limits$max_vertices`/`max_parts`. `validate_coordinate_bounds()` checks the bbox against `bounds_check$xmin/xmax/ymin/ymax`, treating any violation as an error (unlike `validate_wkt_geometry()`'s inline bounds check, which is warning-only). `validate_geographic_context()` additionally checks for physically-impossible lon/lat values (`|lon| > 180`, `|lat| > 90`) in `strict_mode`, and warns if the bbox collapses to a single point. `validate_projected_context()` warns if a "projected" geometry is actually lon/lat, or if coordinate magnitudes exceed `1e8` (possible wrong-CRS indicator). `get_validation_defaults()` returns a `switch()`-selected bundle of `area_limits`/`complexity_limits`/`bounds_check` for `"geographic"`, `"projected"`, `"custom"`, or a final fallback context. `calculate_complexity_score()` is a `log10(vertices) + log10(parts)`, scaled by a `geometry_type`-keyed multiplier (`POINT` 1.0 up to `GEOMETRYCOLLECTION` 2.5), rounded to 2 decimals - it is the only one of this family actually called by another function in the file (`validate_geometry_complexity()`). None of these eight functions is called by `validate_wkt_geometry()` or by any other file in the package.

#### Property-name / synonym resolution

##### `validate_properties(properties, property_lookup = "ssurgo", strict_mode = TRUE, max_invalid_pct = 50, performance_threshold = 10)`

**Purpose**: Simple, source-agnostic validation of a vector of property names against a lookup of "known" names, without synonym matching.

**Returns**: List with `valid`, `errors`, `warnings`, `property_stats` (`total_requested`, `valid_properties`, `invalid_properties`, `validation_method`, `invalid_property_names`).

**Algorithm/behavior**: Rejects non-character/empty `properties` immediately. Resolves the available-property set via `get_available_properties(property_lookup)`; if that returns zero properties, validation is skipped entirely (a warning is recorded, not an error). Otherwise computes `valid_props <- intersect(properties, available_props)` and `invalid_props <- setdiff(...)`; if `strict_mode` and the invalid fraction exceeds `max_invalid_pct`, marks `valid <- FALSE`. Separately warns if the *valid* property count exceeds `performance_threshold`, and fails if there are zero valid properties regardless of `strict_mode`. This is the function actually called elsewhere in the package (`multivariate-adjustment.R`, `validation-diagnostics.R`) - always with `property_lookup = "laboratory"` in current call sites.

##### `get_available_properties(property_lookup)`

**Purpose**: Normalizes four different kinds of `property_lookup` input into a plain character vector of available property names.

**Algorithm/behavior**: Dispatches on the type of `property_lookup`: a length-1 character string is treated as a named source and passed to `get_predefined_properties()`; a function is called and its result coerced via `as.character()`; a longer character vector is returned as-is; a data frame is passed to `extract_properties_from_dataframe()`; anything else logs a warning and returns `character(0)`. The whole dispatch is wrapped in a `tryCatch` that degrades to `character(0)` with a warning on any error.

##### `get_predefined_properties(source_name)`

**Purpose**: Named canned property lists for `"ssurgo"`, `"nrcs"`, `"laboratory"`, and `"basic"` sources.

**Algorithm/behavior**: `"nrcs"`, `"laboratory"`, and `"basic"` return small hard-coded character vectors. `"ssurgo"` is different in kind: it calls `create_ssurgo_property_lookup_working()` - a function defined not in `utils.R` but in `R/ssurgo-acquisition.R` - inside a `tryCatch`, and extracts its `$Property` column if present. This is the one place `utils.R` reaches outside its own file for a live function call (see Known Limitations); the call is defensively guarded, so if `ssurgo-acquisition.R`'s function isn't loaded/available the branch just logs a WARN and returns `character(0)` rather than erroring. Unrecognized `source_name` also logs a WARN and returns `character(0)`.

##### `extract_properties_from_dataframe(df)`

**Purpose**: Pulls a property-name vector out of a data frame by checking a fixed list of likely column names (`property`, `Property`, `name`, `variable`, `column`) in order, falling back to the first column if none match (with a DEBUG log noting the fallback).

##### `validate_properties_with_synonyms(properties, property_lookup = "ssurgo", strict_mode = TRUE, max_invalid_pct = 50, performance_threshold = 10)`

**Purpose**: The enhanced, synonym-aware property validator - the one actually used by SSURGO acquisition, Monte Carlo, and statistics code paths - that resolves names like `"clay"`/`"Clay%"`/`"claytotal"` against a fixed internal SSURGO synonym table.

**Returns**: List with `valid`, `errors`, `warnings`, `property_stats` (`total_requested`, `valid_properties`, `invalid_properties`, `success_rate`, `synonym_matches`, `validation_method`, `valid_property_names`, `invalid_property_names`).

**Algorithm/behavior**: Unlike `validate_properties()`, this function ignores its `property_lookup` argument for matching purposes (it is accepted but not consulted) and instead matches every input name against a **hard-coded** `ssurgo_properties` list of ~13 canonical SSURGO property names, each mapped to a small vector of synonyms/variants (e.g. `"claytotal" = c("claytotal", "clay_total", "clay", "claytotal_r")`). For each requested property, it first checks direct membership in any synonym vector; failing that, it strips a trailing `_r`/`_l`/`_h` suffix and retries the same membership check against the stripped name. Unmatched names accumulate as `invalid_properties`. `invalid_percentage > max_invalid_pct` sets `valid <- FALSE` with an error; a nonzero-but-acceptable invalid count instead produces a warning listing the invalid names. A separate warning fires if `total_properties > performance_threshold`. Because the synonym table here is separate from (and structured oppositely to) `get_default_synonyms()`'s table, the two can drift independently.

##### `get_default_synonyms(property_lookup)`

**Purpose**: Returns a named-list synonym table (canonical name -> vector of synonyms) for `"ssurgo"`, `"nrcs"`, or `"laboratory"`, or an empty list for anything else. Structurally similar in spirit to, but not the same data as, the table hard-coded inside `validate_properties_with_synonyms()`, and not called by it or by any other function in the file.

##### `create_property_lookup(properties, synonyms = NULL, metadata = NULL)`

**Purpose**: Builds a small `"property_lookup"`-classed object bundling a property vector, optional synonym table, and metadata, with two closures attached: `$get_properties()` (returns `lookup$properties`) and `$validate(props, strict_mode = TRUE)`.

**Algorithm/behavior**: `$validate()` dispatches to `validate_properties(props, lookup$properties, strict_mode)` when `lookup$synonyms` is `NULL`, or to `validate_properties_with_synonyms(props, lookup$properties, strict_mode = strict_mode)` otherwise. A source comment at this call site notes a previously-existing bug: the synonym-matching branch used to pass `lookup$synonyms` and `strict_mode` positionally, which landed in `validate_properties_with_synonyms()`'s `strict_mode` and `max_invalid_pct` argument slots respectively (since that function has no `synonyms` parameter at all) - silently mis-validating whenever a lookup object had synonyms set. The current code passes `strict_mode` by name to avoid this.

### Configuration Management

#### `get_default_configuration(config_type = "full")`

**Purpose**: Returns the package-wide baseline configuration.

**Parameters**: `config_type` - `"full"` (everything), `"minimal"` (just `data_processing`, `monte_carlo$n_simulations`, `validation$run_validation`), or `"validation"` (just the `validation` sub-list).

**Returns**: A nested list with top-level sections `workflow`, `data_processing` (`max_depth = 250`, `exclude_unsuitable_horizons = TRUE`, `depth_units = "cm"`, `missing_value_strategy = "interpolate"`), `monte_carlo` (`n_simulations = 100`, `distribution_type = "triangular"`, `confidence_level = 0.95`, `random_seed = 12345`), `correlation` (`method = "pearson"`, `minimum_samples = 30`, `handle_missing = "pairwise.complete.obs"`), `gp_models` (`grouping_strategy = "auto"`, `min_profiles_per_group = 3`, `min_observations_per_group = 15`, `optimize_hyperparameters = TRUE`), `validation` (`run_validation = TRUE`, `quality_threshold = 0.7`, `generate_plots = TRUE`, `max_outlier_percentage = 5`), and `logging` (`log_level = "INFO"`, `log_file = NULL`, `include_timestamp = TRUE`).

#### `load_configuration(config_path, config_type = "auto", validate_config = TRUE)`

**Purpose**: Loads a JSON or R-source configuration file from disk, merges it over the package defaults, and optionally validates it.

**Algorithm/behavior**: If `config_path` doesn't exist, logs a WARN and returns `get_default_configuration()` untouched (no error). Otherwise auto-detects `config_type` from the file extension (`json`/`r`, defaulting to `json` for anything unrecognized) unless explicitly given. Loads via `jsonlite::fromJSON(config_path, simplifyVector = FALSE)` for JSON, or `source(config_path)$value` for R (i.e. an R script whose last top-level expression is treated as the config's `value`); any other `config_type` is a hard `stop()`. The loaded config is deep-merged onto `get_default_configuration()` via `merge_configurations()`. If `validate_config`, runs the internal (unexported) `validate_configuration()` and logs a WARN (not an error) listing any validation errors.

#### `merge_configurations(default_config, user_config, deep_merge = TRUE)`

**Purpose**: Layers a user-supplied configuration over the package defaults.

**Algorithm/behavior**: Returns `default_config` unchanged if `user_config` is `NULL`. With `deep_merge = FALSE`, does a shallow, top-level-only overwrite (`merged_config[[name]] <- user_config[[name]]` for each name in `user_config`). With `deep_merge = TRUE` (the default), uses a local recursive helper `merge_recursive()`: for any name present in both `default` and `user` where both sides are themselves lists, it recurses; otherwise the user's value simply replaces the default's, whole. Non-list leaf values (e.g. a user override of a single number) always take the user value directly, even under deep merge.

#### `validate_parameters(parameters, parameter_specs, strict_mode = TRUE)`

**Purpose**: Generic parameter validation against a spec list - required-ness, type, numeric range, allowed choices, and length - used ahead of expensive/side-effecting operations (e.g. SSURGO downloads).

**Returns**: List with `valid`, `errors`, `warnings`, `corrected_parameters` (currently always identical to the input `parameters` - nothing in this function actually mutates or clips values despite the name).

**Algorithm/behavior**: Iterates `names(parameter_specs)`. A missing/zero-length parameter is only an error if `param_spec$required` is `TRUE` (via `%||% FALSE`); otherwise the spec is skipped entirely for that parameter (`next`), meaning type/range/choices/length checks never run against absent optional parameters. Type checking (`param_spec$type` in `c("character","numeric","integer","logical","list")`) uses `is.*()` predicates, with `"integer"` accepting either a true integer or a numeric whose values are all whole numbers; an unrecognized `type` string is treated as automatically valid. Range checking (`param_spec$range`, a 2-element `c(min,max)`) applies only to numeric values and is vectorized (`any(out_of_range)`); in `strict_mode` an out-of-range value is an error, otherwise a warning. Choices checking (`param_spec$choices`) branches on whether `param_value` has length 1 (single value must be `%in%` choices) or is a vector (each element checked, invalid ones collected and reported together) - both branches always produce an error on any mismatch, irrespective of `strict_mode`. `min_length`/`max_length` checks are independent of the choices/range/type checks; `max_length` violations are errors in `strict_mode` and warnings otherwise, while `min_length` violations are always errors.

#### `export_workflow_metadata(workflow_results, output_path, include_data_summary = TRUE)`

**Purpose**: Serializes a JSON metadata snapshot of a completed workflow run (versions, execution summary, configuration, quality assessment, optional data summary) to disk.

**Algorithm/behavior**: Assembles `workflow_info` (timestamp, `R.version.string`, platform, and a filtered `installed.packages()` lookup for `dplyr`/`tidyr`/`terra`/`soilDB`/`readr`/`jsonlite` via the internal `get_package_versions()`), plus `execution_summary`, `configuration`, and `quality_assessment` sections built by internal `extract_*()` helpers that defensively look for specific named elements (`processing_metadata`, `validation_results`, `quality_report`, `processed_data`) inside `workflow_results` and fall back to `"Not available"` placeholder strings wherever the expected shape isn't found. If `include_data_summary`, also attaches a `data_summary` (row/column counts and up to 5 `_r`-suffixed numeric property summaries: mean and missing-percentage). Writes the assembled list via `jsonlite::write_json(metadata, output_path, pretty = TRUE, auto_unbox = TRUE)` inside a `tryCatch`, returning `TRUE`/`FALSE` for success/failure (logged either way). Not currently called elsewhere in the package (see Known Limitations).

### Logging & Error Handling

#### `setup_logging(log_file = NULL, log_level = "INFO", include_timestamp = TRUE, max_log_size = 10)`

**Purpose**: Initializes the package-wide logging configuration, stashed in `options()` for `log_message()` to read on every call.

**Algorithm/behavior**: Creates the log file's parent directory (if `log_file` is non-`NULL`) via `dir.create(..., recursive = TRUE, showWarnings = FALSE)`. Builds a `log_config` list (`log_file`, `log_level`, `include_timestamp`, `max_log_size`, `start_time = Sys.time()`) and stores it with `options(soil_workflow_log_config = log_config)` - a single global, process-wide logging configuration (not scoped per-call or per-package-instance). Logs its own initialization via `log_message()`, then returns `log_config` invisibly... (actually visibly; the return itself is not wrapped in `invisible()`).

#### `log_message(level, message, category = NULL, details = NULL)`

**Purpose**: The central logging call used throughout the entire package - console + optional rotating log file, level-filtered.

**Algorithm/behavior**: Reads the current `soil_workflow_log_config` option (falling back to an ad hoc `INFO`/console-only default if `setup_logging()` was never called). Compares `level` against `log_config$log_level` using a fixed `c(DEBUG=1, INFO=2, WARN=3, ERROR=4)` hierarchy; messages below the configured threshold (or with an unrecognized `level` string, since `NA` fails the numeric comparison) are silently dropped (`invisible()` early return - an unrecognized level is dropped rather than erroring or defaulting to always-shown). Builds the formatted line as `[timestamp] [LEVEL] message [category] - details`, with the timestamp and the `category`/`details` suffixes each optional. `WARN`/`ERROR` go to `stderr()`; everything else goes to stdout via plain `cat()`. If a log file is configured, checks its current size against `max_log_size` (MB) and calls the internal `rotate_log_file()` (renaming the existing file with a timestamp suffix) before appending the new line.

#### `handle_workflow_error(error, context = "Unknown", recovery_action = "stop", max_retries = 3)`

**Purpose**: Centralized error handling with a configurable recovery strategy, used inside essentially every `tryCatch(..., error = function(e) ...)` block across the package.

**Algorithm/behavior**: Always logs an ERROR line (`"Error in <context> : <error$message>"`) and, if `error$call` is present, a DEBUG line with the deparsed call. Then a `switch()` on `recovery_action`: `"stop"` logs and re-raises via `stop(error_message, call. = FALSE)` (the original condition object/class is not preserved - callers see a fresh error rather than the original condition); `"warn"` logs, raises an R `warning()`, and *also* returns `list(status = "warning", message = ...)` (so callers using this mode must not expect the function to halt - it both warns and returns a value); `"continue"` just returns `list(status = "error_ignored", message = ...)`; `"retry"` returns `list(status = "retry", message = ..., max_retries = max_retries)` - note this function only ever reports a retry request, it does not itself loop or re-invoke anything; any actual retry loop must live in the caller. An unrecognized `recovery_action` string falls through to the `switch()`'s default branch, which logs and `stop()`s exactly like `"stop"`.

#### `track_progress(current, total, message = "Processing", update_frequency = 10)`

**Purpose**: Throttled progress logging plus a rough time-remaining estimate for long-running loops (downloads, Monte Carlo simulations).

**Algorithm/behavior**: Skips logging entirely unless `current` is a multiple of `update_frequency` or equals `total` (so the very first calls, other than `current == 1` for timer setup, are typically not logged unless `update_frequency` is 1). Computes `percentage <- round(current/total * 100, 1)`. If `current > 0`, reads a `progress_start_time` option (defaulting to `Sys.time()` - i.e. effectively zero elapsed time - if never set), computes `rate <- current/elapsed` and `remaining <- (total - current)/rate`, and appends a `"X min remaining"`/`"X sec remaining"` string (minutes once remaining exceeds 60 seconds). On `current == 1`, sets `options(progress_start_time = Sys.time())` for subsequent calls to read - but this reset happens *after* the current call's own elapsed-time calculation above, and only when `current == 1` passes the `update_frequency` gate at the top of the function (e.g. it will be skipped entirely if `update_frequency` isn't 1 or a divisor that happens to allow row 1 through, unless `total == 1` too). Logs via `log_message("INFO", ..., category = "Progress")`.

### General Data/Stats Helpers

#### `safe_coalesce(data, column_names, target_type = "numeric", default_value = NA)`

**Purpose**: Null-safe coalescing of a value across several candidate columns, in priority order, with type coercion.

**Algorithm/behavior**: Restricts to columns that actually exist in `data` (logging a WARN and returning an all-`default_value` vector if none do). Starting from a `result` vector pre-filled with `default_value`, walks the existing columns in the caller-supplied priority order; each column is coerced to `target_type` (`numeric` via `suppressWarnings(as.numeric())`, `character` via `as.character()`, `logical` via `as.logical()`) and its values are written into any positions of `result` still needing a value - `is.na(result)` for numeric targets, or `is.na(result) | result == "" | result == default_value` for the other two target types. Logs the count of values ultimately filled.

#### `handle_missing_values(data, strategy = "interpolate", columns = NULL, group_by = NULL)`

**Purpose**: Applies one of several missing-value strategies across a set of columns, optionally per-group.

**Parameters**: `strategy` - `"remove"`, `"interpolate"` (default), `"mean"`, `"median"`, or `"forward_fill"` (see the internal vector-level helper below). `columns` - defaults to all numeric columns in `data` if `NULL`. `group_by` - optional grouping column names for per-group imputation.

**Algorithm/behavior**: For each target column with at least one `NA`, if `group_by` is supplied and all its names exist in `data`, applies the strategy per-group via a `dplyr` pipe (`group_by() |> mutate(!!col := handle_missing_values_vector(.data[[col]], strategy)) |> ungroup()`); otherwise applies the internal `handle_missing_values_vector()` helper globally on the whole column. Logs, per column, how many of the originally-missing values were filled.

#### `detect_outliers(data, method = "iqr", threshold = NULL, return_indices = FALSE)`

**Purpose**: Dispatching entry point for outlier detection across three methods.

**Algorithm/behavior**: If `threshold` is `NULL`, picks a method-specific default (`1.5` for `"iqr"`, `3` for `"zscore"`, `3.5` for `"modified_zscore"`, `1.5` as the catch-all). Dispatches via `switch()` to one of three internal (unexported) helpers - `detect_outliers_iqr()` (`x < Q1 - k*IQR | x > Q3 + k*IQR`), `detect_outliers_zscore()` (`abs((x-mean)/sd) > k`), or `detect_outliers_modified_zscore()` (`abs(0.6745*(x-median)/mad) > k`) - each returning a logical vector the same length as the input. If `return_indices`, converts the logical vector to `which()` indices before returning.

#### `normalize_values(x, method = "minmax", center = TRUE, scale = TRUE)`

**Purpose**: Rescales a numeric vector by one of four methods.

**Algorithm/behavior**: Returns `x` unchanged if every value is `NA`. `"minmax"` maps to `[0,1]`, special-casing a constant vector (`min == max`) to a vector of `0.5`s rather than dividing by zero. `"zscore"` uses base `scale(x, center, scale)` when both `center` and `scale` are `TRUE`, or manual centering/scaling only when just one is requested, or returns `x` unchanged if neither. `"robust"` uses median/MAD instead of mean/SD, again special-casing `mad == 0` to avoid division by zero (returns `x - median` in that case). `"quantile"` computes rank-based normalization: `(rank(x, na.last = "keep") - 1) / (n_valid - 1)`. Any other `method` string is a hard `stop()`.

#### `safe_correlation(x, y = NULL, method = "pearson", handle_constant = "warn")`

**Purpose**: Numerically defensive correlation computation, either pairwise (`x`, `y` both vectors) or as a full correlation matrix (`x` a data frame/matrix, `y = NULL`).

**Algorithm/behavior - matrix mode**: Requires `x` be a data frame or matrix (`stop()`s otherwise); if a data frame, drops non-numeric columns first. Flags constant columns (`var(col, na.rm=TRUE) == 0` or all-`NA`) and either drops them (`handle_constant = "remove"`) or just warns (`"warn"`) while keeping them. If fewer than 2 columns remain, returns a trivial `1x1` matrix of `1` rather than erroring. Computes `cor(x, method, use = "pairwise.complete.obs")`; if the result has any `NA` and `method == "pearson"`, automatically retries with `method = "spearman"` before returning. Any error inside the `cor()` call is caught and replaced with an identity matrix (`diag(ncol(x))`) sized to match, rather than propagating.

**Algorithm/behavior - pairwise mode**: `stop()`s if `x`/`y` differ in length. Restricts to `complete.cases(x, y)`; if fewer than 3 complete pairs remain, warns and returns `NA`. If either variable is constant (`var == 0`), returns `0` for `handle_constant %in% c("zero","warn")` (logging in the `"warn"` case) or `NA` otherwise. Otherwise computes `cor(x_clean, y_clean, method)` inside a `tryCatch`, returning `NA` on any computation error.

#### `calculate_confidence_intervals(data, statistic = "mean", confidence_level = 0.95, method = "normal", n_bootstrap = 1000)`

**Purpose**: Confidence interval computation for a mean, median, or proportion, by normal-approximation, t-distribution, or bootstrap.

**Parameters**: `statistic` - `"mean"`, `"median"`, or `"proportion"` (proportion is treated identically to mean - `mean(x)` - throughout; there is no dedicated binomial/Wilson-type interval logic). `method` - `"normal"` (z-based), `"t"` (t-based, uses `df = n-1`), or `"bootstrap"` (percentile bootstrap over `n_bootstrap` resamples).

**Algorithm/behavior**: Drops `NA`s from `data` up front (logging the count removed - note the logged count is computed *after* the removal, from the now-`NA`-free vector, so it will always log `0` removed regardless of how many were actually stripped). Returns an all-`NA` result immediately if nothing is left. For `"bootstrap"`, draws `n_bootstrap` resamples of `sample(data, length(data), replace = TRUE)`, computes the requested `statistic` on each resample via `replicate()`, and takes the `alpha/2`/`1-alpha/2` sample quantiles of the bootstrap distribution as the CI bounds (an unrecognized `statistic` inside the bootstrap loop is a hard `stop()`). For `"normal"`/`"t"`, computes `se <- sd(data)/sqrt(n)` and a margin of error from either `qnorm(1 - alpha/2)` or `qt(1 - alpha/2, df = n-1)`, applied symmetrically around the sample mean (note: for `"t"`/`"normal"`, `estimate` is always the mean regardless of the requested `statistic` - `"median"`/`"proportion"` are only actually honored in the `"bootstrap"` branch). Any other `method` value is a hard `stop()`. Returns `lower`, `upper`, `estimate`, `confidence_level`, `method`.

### File I/O

#### `read_soil_data(file_path, file_type = "auto", validate_on_load = TRUE, sheet_name = NULL)`

**Purpose**: Format-agnostic soil data reader with optional post-load quality validation.

**Algorithm/behavior**: `stop()`s if `file_path` doesn't exist. Auto-detects type from the file extension (`tools::file_ext()`) when `file_type = "auto"`, defaulting to `"csv"` for any unrecognized extension. Dispatches to the internal helpers `read_csv_robust()` (tries `readr::read_csv()`, falls back to base `read.csv()` on error), `read_excel_robust()` (requires `readxl`; reads a named sheet or the default), `readRDS()` directly, or `read_tsv_robust()` (tries `readr::read_tsv()`, falls back to `read.delim()`); any other `file_type` is a hard `stop()`. If `validate_on_load`, runs `validate_data_quality(data)` (with no `required_columns`/`numeric_columns` supplied, so only the missing-data and consistency sub-checks are meaningful) and logs a WARN if the resulting score is below 0.7 - this is advisory only, it never blocks the load or alters the returned data.

#### `write_soil_data(data, file_path, file_type = "csv", include_metadata = TRUE, create_backup = TRUE)`

**Purpose**: Format-agnostic soil data writer with an optional pre-write backup of any existing file at the destination.

**Algorithm/behavior**: If `create_backup` and a file already exists at `file_path`, copies it aside first via the internal `create_backup_filename()` (timestamped, same directory, same extension) - note this duplicates, with an independent implementation, the timestamped-filename logic also present inline in `backup_data()` (see Known Limitations). Ensures the destination directory exists (`dir.create(..., recursive = TRUE, showWarnings = FALSE)`). Writes via `readr::write_csv()`, the internal `write_excel_with_metadata()` (requires `openxlsx`; writes a `"Data"` sheet and, if `include_metadata`, a second `"Metadata"` sheet with export date/row count/column count/R version), or `saveRDS()`; any other `file_type` is a hard `stop()`. The whole write is wrapped in `tryCatch`, returning `TRUE` on success or `FALSE` (after logging an ERROR) on failure - failures are swallowed into a boolean return rather than propagated as R conditions.

#### `backup_data(source_path, backup_dir = NULL, max_backups = 5)`

**Purpose**: Creates a timestamped backup copy of a file with automatic pruning of older backups.

**Algorithm/behavior**: `stop()`s if `source_path` doesn't exist. Defaults `backup_dir` to the source file's own directory. Builds a timestamped filename (`<base>_backup_<YYYYMMDD_HHMMSS>.<ext>`) **inline**, independently of the very similar internal `create_backup_filename()` helper used by `write_soil_data()` (a small code-duplication - see Known Limitations). Copies the file via `file.copy()`; on success, logs and calls the internal `cleanup_old_backups()` (which lists files matching `<base>_backup_` in `backup_dir`, sorts by modification time descending, and `file.remove()`s everything beyond the newest `max_backups`), returning the new backup's path. On copy failure, `stop("Failed to create backup")`.

## Internal Connections

Calls verified directly against the source (only functions defined in this file; `log_message()` calls are omitted from most nodes below for brevity since nearly every function calls it):

```
utils.R
│
├── is_unsuitable()                         [no internal calls beyond log_message()]
├── validate_data_quality()                 [does NOT call detect_outliers() - inlines its own IQR check]
├── check_required_columns()                [no internal calls; unused elsewhere]
├── validate_numeric_ranges()               [no internal calls]
│
├── standardize_property_names()
│   └── get_default_property_mapping()      (internal helper)
├── convert_depth_units()                   [no internal calls; unused elsewhere]
│
├── safe_coalesce()                         [no internal calls]
├── handle_missing_values()
│   └── handle_missing_values_vector()      (internal helper)
│
├── read_soil_data()
│   ├── read_csv_robust() / read_excel_robust() / read_tsv_robust()   (internal helpers)
│   └── validate_data_quality()
├── write_soil_data()
│   ├── create_backup_filename()            (internal helper)
│   └── write_excel_with_metadata()         (internal helper)
├── backup_data()
│   ├── (inline timestamped-filename logic - does NOT call create_backup_filename())
│   └── cleanup_old_backups()               (internal helper)
│
├── load_configuration()
│   ├── get_default_configuration()
│   ├── merge_configurations()
│   └── validate_configuration()            (internal, unexported)
├── get_default_configuration()             [no internal calls]
├── merge_configurations()                  [no internal calls; self-contained recursive helper]
├── validate_parameters()                   [no internal calls]
├── export_workflow_metadata()
│   ├── get_package_versions()              (internal helper)
│   ├── extract_execution_summary()         (internal helper)
│   ├── extract_configuration_info()        (internal helper)
│   ├── extract_quality_summary()           (internal helper)
│   └── extract_data_summary()              (internal helper)
│
├── safe_correlation()                      [no internal calls]
├── calculate_confidence_intervals()        [no internal calls]
├── normalize_values()                      [no internal calls]
├── detect_outliers()
│   └── detect_outliers_iqr() / detect_outliers_zscore() / detect_outliers_modified_zscore()  (internal helpers)
│
├── setup_logging()                         [calls log_message() only]
├── log_message()
│   └── rotate_log_file()                   (internal helper; only if log file exceeds max_log_size)
├── handle_workflow_error()                 [calls log_message() only]
├── track_progress()                        [calls log_message() only]
│
├── validate_properties()
│   └── get_available_properties()
│       ├── get_predefined_properties()
│       │   └── create_ssurgo_property_lookup_working()   *** defined in R/ssurgo-acquisition.R, NOT utils.R ***
│       └── extract_properties_from_dataframe()
├── validate_properties_with_synonyms()     [own hard-coded synonym table; does NOT call validate_properties()
│                                             or get_available_properties()]
├── get_default_synonyms()                  [no internal calls; separate synonym table from the one above]
├── create_property_lookup()                [closure dispatches to validate_properties() or
│                                             validate_properties_with_synonyms() at call time]
│
├── validate_wkt_geometry()                 [inlines all parsing/bbox/area/bounds logic itself;
│                                             does NOT call any of the granular functions below]
├── validate_wkt_string()                   [no internal calls; unused elsewhere]
├── parse_wkt_geometry()                    [no internal calls other than log_message(); unused elsewhere]
├── validate_geometry_validity()            [no internal calls; unused elsewhere]
├── validate_geometry_area()                [no internal calls; unused elsewhere]
├── validate_geometry_complexity()
│   └── calculate_complexity_score()        (the only cross-call within this granular geometry family)
├── validate_coordinate_bounds()            [no internal calls; unused elsewhere]
├── get_validation_defaults()               [no internal calls; unused elsewhere]
├── validate_geographic_context()           [no internal calls; unused elsewhere]
└── validate_projected_context()            [no internal calls; unused elsewhere]
```

## Dependencies

### External packages
- `dplyr` - `group_by()`/`mutate()`/`ungroup()` (grouped missing-value handling in `handle_missing_values()`); `across()`/`all_of()`.
- `tidyr` - imported per `DESCRIPTION`; not directly called inside `utils.R` itself (its use lives in other package files).
- `readr` - `read_csv()`, `read_tsv()`, `write_csv()` (with base-R `read.csv()`/`read.delim()` fallbacks on error).
- `jsonlite` - `fromJSON()` (config loading), `write_json()` (metadata export).
- `tools` - `file_ext()`, `file_path_sans_ext()` (file-type/backup-filename detection throughout File I/O).
- `terra` - `vect()`, `is.valid()`, `crds()`, `expanse()`, `geomtype()`, `ext()`, `nrow()`, `geom()`, `is.lonlat()`, `makeValid()`, `crs()` - used only inside the WKT/geometry-validation family; `requireNamespace("terra", ...)`-guarded in `validate_wkt_geometry()`.
- `readxl` (Suggests) - guarded with `requireNamespace()` inside `read_excel_robust()`; reading Excel files errors with a clear message if absent.
- `openxlsx` (Suggests) - guarded with `requireNamespace()` inside `write_excel_with_metadata()`.
- Base R - `stats::sd/var/cor/quantile/qnorm/qt/median/mad/rank`, `installed.packages()`, file/dir functions (`file.exists`, `file.copy`, `file.remove`, `dir.create`, `file.info`, `file.size`, `file.rename`), `Sys.time()`, `options()`/`getOption()`.

### soilSIM-internal dependencies
`utils.R` sources nothing else in the package and is intended as a pure leaf module - with one verified exception: `get_predefined_properties("ssurgo")` calls `create_ssurgo_property_lookup_working()`, which is defined in `R/ssurgo-acquisition.R`. This call is defensively wrapped in `tryCatch()` (degrading to a `WARN` + empty vector if unavailable), so it does not create a hard circular *load-order* dependency, but it does mean `utils.R` is not a fully self-contained leaf at the call-graph level. See Known Limitations.

### Consumers (functions from this file used elsewhere, per source grep)
- **`is_unsuitable()`** - the single most widely reused function: `R/ssurgo-acquisition.R`, `R/ssurgo-processing.R`, `R/data-infilling.R` (three call sites), `R/gp-modeling.R`, `R/monte-carlo.R`.
- **`validate_data_quality()`** - `R/ssurgo-acquisition.R`, `R/ssurgo-processing.R`, `R/gp-modeling.R`, `R/multivariate-adjustment.R`, `R/statistics.R`, `R/validation-diagnostics.R` (multiple call sites).
- **`validate_properties_with_synonyms()`** - `R/ssurgo-acquisition.R`, `R/monte-carlo.R`, `R/statistics.R`, `R/data-infilling.R`.
- **`validate_wkt_geometry()`** - `R/ssurgo-acquisition.R`.
- **`validate_properties()`** (the non-synonym variant) - `R/multivariate-adjustment.R`, `R/validation-diagnostics.R` (all current call sites pass `"laboratory"` as the lookup).
- **`validate_numeric_ranges()`** - `R/ssurgo-processing.R`, `R/validation-diagnostics.R`.
- **`standardize_property_names()`** - `R/gp-modeling.R`, `R/ssurgo-processing.R` (two call sites).
- **`safe_coalesce()`** - `R/gp-modeling.R`.
- **`handle_missing_values()`** - `R/statistics.R`.
- **`detect_outliers()`** - `R/monte-carlo.R` (two call sites), `R/statistics.R` (two call sites), `R/ssurgo-processing.R`, `R/validation-diagnostics.R` (three call sites).
- **`safe_correlation()`** - `R/statistics.R` (multiple call sites), `R/validation-diagnostics.R` (three call sites).
- **`calculate_confidence_intervals()`** - `R/statistics.R`.
- **`get_default_configuration()`** - `R/gp-modeling.R`, `R/multivariate-adjustment.R` (multiple call sites), `R/ssurgo-acquisition.R` (two call sites), `R/statistics.R`, `R/validation-diagnostics.R` (multiple call sites), `R/monte-carlo.R`.
- **`merge_configurations()`** - `R/monte-carlo.R` (two call sites, one inside its own roxygen-documented helper), `R/statistics.R` (two call sites), `R/ssurgo-processing.R`.
- **`load_configuration()`** - `R/ssurgo-acquisition.R`.
- **`setup_logging()`** - `R/monte-carlo.R`, `R/statistics.R`, `R/ssurgo-acquisition.R`.
- **`handle_workflow_error()`** - used pervasively as the `error =` handler inside `tryCatch()` throughout `R/gp-modeling.R`, `R/multivariate-adjustment.R`, `R/monte-carlo.R`, `R/ssurgo-acquisition.R`, `R/validation-diagnostics.R`.
- **`track_progress()`** - `R/multivariate-adjustment.R` (two call sites), `R/monte-carlo.R`, `R/validation-diagnostics.R` (multiple call sites).

Functions defined in this file but **not currently called from any other file in the package** (confirmed by source grep): `check_required_columns()`, `convert_depth_units()`, `export_workflow_metadata()`, `validate_wkt_string()`, `parse_wkt_geometry()`, `validate_geometry_validity()`, `validate_geometry_area()`, `validate_geometry_complexity()`, `validate_coordinate_bounds()`, `get_validation_defaults()`, `validate_geographic_context()`, `validate_projected_context()`, `get_default_synonyms()`, `get_available_properties()`/`get_predefined_properties()`/`extract_properties_from_dataframe()`/`create_property_lookup()` (these four are called only from within `utils.R` itself, via `validate_properties()`).

## Data Flow In/Out

`utils.R` has no SSURGO- or workflow-specific data shape of its own; its functions take generic inputs - data frames, character vectors, WKT strings, numeric vectors, and configuration lists - and return one of a few generic output shapes depending on which function is called:
- **Validation-result lists** with a consistent `valid`/`errors`/`warnings` shape (`validate_data_quality()`, `validate_properties()`, `validate_properties_with_synonyms()`, `validate_wkt_geometry()`, `validate_wkt_string()`, `validate_parameters()`, `check_required_columns()`, `validate_numeric_ranges()`, and the granular geometry validators), often with an additional `*_stats`/`property_stats` sub-list of diagnostic numbers.
- **Logical vectors** (`is_unsuitable()`, `detect_outliers()` when `return_indices = FALSE`) or **integer index vectors** (`detect_outliers()` when `return_indices = TRUE`).
- **Cleaned/transformed data frames** (`handle_missing_values()`, `standardize_property_names()`, `validate_numeric_ranges()$corrected_data`).
- **Numeric vectors/matrices** (`normalize_values()`, `convert_depth_units()`, `safe_coalesce()`, `safe_correlation()`).
- **Small structured lists** for statistics (`calculate_confidence_intervals()`) and configuration (`get_default_configuration()`, `load_configuration()`, `merge_configurations()`).
- **Side effects only** (`setup_logging()`, `log_message()`, `track_progress()` write to console/file/`options()` and return `invisible()` or the config list itself; `write_soil_data()`/`backup_data()`/`export_workflow_metadata()` write to disk and return a path or a success boolean).

## Known Limitations

- **`get_predefined_properties("ssurgo")` breaks the leaf-module property**: it calls `create_ssurgo_property_lookup_working()`, which is defined in `R/ssurgo-acquisition.R`, not `utils.R`. The call is `tryCatch`-guarded (degrades to `WARN` + `character(0)` if unavailable), so it is not a hard failure, but it means `utils.R` is not a fully self-contained leaf at the call-graph level despite otherwise sourcing nothing else in the package.
- **Two independent, non-communicating WKT validation code paths**: `validate_wkt_geometry()` (the one actually wired into `R/ssurgo-acquisition.R`) inlines all of its own parsing/bbox/area/bounds logic and does not call `parse_wkt_geometry()`, `validate_geometry_validity()`, `validate_geometry_area()`, `validate_geometry_complexity()`, `validate_coordinate_bounds()`, `get_validation_defaults()`, `validate_geographic_context()`, or `validate_projected_context()` - all nine of which are exported but otherwise dead code (unused anywhere in the package). `validate_wkt_string()` is a third, simpler string-only checker, also unused elsewhere.
- **Two independent SSURGO synonym tables that can drift**: `validate_properties_with_synonyms()` hard-codes its own ~13-entry SSURGO synonym table inline, structurally unrelated to (and not sharing data with) `get_default_synonyms("ssurgo")`'s separate table or `get_default_property_mapping("ssurgo")`'s (used by `standardize_property_names()`). All three tables encode overlapping but independently-maintained SSURGO name/synonym knowledge.
- **Duplicated backup-filename logic**: `write_soil_data()` uses the internal helper `create_backup_filename()`, while `backup_data()` builds an equivalent timestamped filename inline with its own independent code rather than calling that helper.
- **`validate_data_quality()`'s duplicate-source-of-truth for the pass/fail threshold**: the function reads `quality_thresholds$min_completeness` for the final `validation_passed` check, but falls back to a literal `0.7` if unset, whereas the auto-generated default `quality_thresholds$min_completeness` (when the caller passes no thresholds at all) is `0.8` - the two default values are inconsistent with each other.
- **`calculate_confidence_intervals()`'s missing-value log message always reports zero removed**: the "removed N missing values" count is computed from `data` *after* it has already had its `NA`s stripped, so the logged count is always `0` regardless of how many values were actually dropped.
- **`export_workflow_metadata()` and `check_required_columns()` are unused elsewhere** in the current package (confirmed by source grep) - both are fully implemented and exported but have no call sites outside their own definitions.
- No `@section Known limitation:` roxygen blocks are present anywhere in `utils.R` itself; the items above were identified by reading the source rather than from author-flagged doc comments (one exception is the `create_property_lookup()` synonym-argument bug, which the file's own source comment documents as a fixed historical bug rather than an open limitation - see that function's write-up above).

## Usage Example

```r
# Validate a set of user-supplied property names against SSURGO synonyms,
# then run a composite data-quality check before doing anything expensive with it.

property_check <- validate_properties_with_synonyms(
  properties = c("clay", "sandtotal_r", "pH", "bogus_property"),
  property_lookup = "ssurgo",
  strict_mode = TRUE
)

if (!property_check$valid) {
  stop("Property validation failed: ", paste(property_check$errors, collapse = "; "))
}

# ssurgo_data is a horizon-level data frame with an `hzname` column
ssurgo_data$unsuitable_horizon <- is_unsuitable(ssurgo_data, hzname_col = "hzname")

quality <- validate_data_quality(
  data = ssurgo_data[!ssurgo_data$unsuitable_horizon, ],
  required_columns = c("cokey", "hzdept_r", "hzdepb_r"),
  numeric_columns = c("claytotal_r", "sandtotal_r", "silttotal_r")
)

log_message("INFO", paste("Data quality grade:", quality$overall_quality$grade), category = "Example")
```
