# Data Acquisition & Processing

## Overview

This functional group is soilSIM’s SSURGO data-source adapter. It is the
entry point of the entire soilSIM pipeline: it acquires raw tabular soil
survey data from USDA-NRCS’s Soil Data Access (SDA) service for a
user-specified area of interest, cleans and standardizes that data
(parsing malformed strings, removing statistical/physical outliers,
aggregating rock-fragment fractions, flagging unsuitable horizons), and
then infills missing low/representative/high (“triplet”) property values
using a hierarchy of pedologically-informed recovery strategies
(horizon-name matching, depth-weighted averaging, within- and
cross-component interpolation, related-property estimation, and
pedotransfer functions such as Saxton-Rawls). The three files covered
here - `ssurgo-acquisition.R`, `ssurgo-processing.R`, and
`data-infilling.R` - together turn a WKT polygon and a list of property
names into a complete, gap-free horizon/component data set with
`_l`/`_r`/`_h` triplet columns, ready to be handed to the Statistics &
Diagnostics module and the Monte Carlo simulation core.

## Core Functions

### 1. **`download_ssurgo_tabular()`** - Master Download Function

**Purpose**: Downloads SSURGO tabular horizon/component/restriction data
for an AOI from Soil Data Access, with input validation, caching,
rock-fragment aggregation, unsuitable-horizon flagging, and data-quality
validation.

**Parameters**:

``` r

download_ssurgo_tabular(
  aoi_wkt,                     # Character. WKT representation of the area of interest.
  properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "ph1to1h2o",
                 "cec7", "om", "wthirdbar", "wfifteenbar"),
  include_restrictions = TRUE, # Logical. Include restriction/horizon-suitability fields.
  cache_dir = NULL,            # Character. Directory for caching; NULL disables caching.
  force_download = FALSE,      # Logical. Bypass cache and force a fresh download.
  validate_data = TRUE,        # Logical. Run data-quality validation on the result.
  verbose = getOption("ssurgo.verbose", FALSE)  # Logical. Emit progress log messages.
)
```

**Returns**: A list with `ssurgo_data` (combined horizon/component data
with restriction flags), `mu` (spatial map-unit data from
[`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)),
`metadata` (download/processing metadata), `validation_results` (from
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
if `validate_data = TRUE`), and `cache_info` (cache hit/write status).

**Algorithm**: Validates inputs via
[`validate_download_inputs_ssurgo()`](https://jjmaynard.github.io/soilSIM/reference/validate_download_inputs_ssurgo.md);
if `cache_dir` is set and not bypassed, checks
[`check_ssurgo_cache()`](https://jjmaynard.github.io/soilSIM/reference/check_ssurgo_cache.md)
and returns cached data immediately on a hit (re-validating if
requested). On a cache miss, it builds the property lookup table
([`create_ssurgo_property_lookup_working()`](https://jjmaynard.github.io/soilSIM/reference/create_ssurgo_property_lookup_working.md)),
resolves the AOI to map unit keys
([`process_aoi_and_get_mukeys_working()`](https://jjmaynard.github.io/soilSIM/reference/process_aoi_and_get_mukeys_working.md)),
executes the SDA SQL query
([`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)),
aggregates rock-fragment volume if `"rfv"` was requested
([`aggregate_rock_fragment_volume_working()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_rock_fragment_volume_working.md)),
adds restriction/unsuitable-horizon indicators
([`add_restriction_indicators_working()`](https://jjmaynard.github.io/soilSIM/reference/add_restriction_indicators_working.md),
or a bare
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)
call when restrictions are excluded), validates the result
([`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)),
writes it to cache if enabled
([`cache_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/cache_ssurgo_data.md)),
and finally assembles timing/quality metadata
([`create_download_metadata()`](https://jjmaynard.github.io/soilSIM/reference/create_download_metadata.md)).

### 2. **`create_ssurgo_property_lookup_working()`**

**Purpose**: Builds the property-name lookup table mapping the 9 default
SSURGO property base names to their `_l`/`_r`/`_h` SSURGO column names.

**Parameters**: None.

**Returns**: A data frame with columns `Property`, `SSURGO_Label_Low`,
`SSURGO_Label_Rep`, `SSURGO_Label_High`, one row per of the 9 known
properties (`sandtotal`, `claytotal`, `silttotal`, `dbovendry`,
`ph1to1h2o`, `cec7`, `om`, `wthirdbar`, `wfifteenbar`).

**Algorithm**: Static
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) construction -
a hardcoded lookup, not derived from the input properties.

### 3. **`process_aoi_and_get_mukeys_working(aoi_wkt, verbose = FALSE)`**

**Purpose**: Converts a WKT area of interest into the set of SSURGO map
unit keys (mukeys) that intersect it.

**Parameters**: `aoi_wkt` - WKT geometry string; `verbose` - emit debug
log messages.

**Returns**: A list with `aoi` (the reprojected `terra` polygon, as a
bounding-box extent polygon), `mu` (the
[`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
raster/spatial result), and `mukey_list` (unique, non-`NA` mukey
integers).

**Algorithm**: `terra::vect(aoi_wkt, crs = "epsg:4326")` →
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)
to `EPSG:5070` (CONUS Albers) → collapses to its bounding-box polygon
via `terra::as.polygons(terra::ext(aoi))` →
`soilDB::mukey.wcs(aoi, db = "gssurgo")` to fetch the gridded map-unit
key coverage → extracts unique non-`NA` mukey values. Both the spatial
conversion and the mukey fetch are wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) blocks that
[`stop()`](https://rdrr.io/r/base/stop.html) with a descriptive message
on failure.

### 4. **`execute_ssurgo_query_working(mukey_list, properties, ssurgo_lookup, include_restrictions = TRUE, verbose = FALSE)`**

**Purpose**: Constructs and executes the SQL query against Soil Data
Access for the requested mukeys and properties.

**Parameters**: `mukey_list` - integer mukeys to filter on;
`properties` - requested property base names; `ssurgo_lookup` - lookup
table from
[`create_ssurgo_property_lookup_working()`](https://jjmaynard.github.io/soilSIM/reference/create_ssurgo_property_lookup_working.md);
`include_restrictions` - whether to join `corestrictions`/`chtexture`
and select restriction fields; `verbose` - log query submission/result
details.

**Returns**: A data frame of raw query results (component, horizon,
requested property `_l/_r/_h` columns, rock-fragment `fragsize_r`, and,
if requested, restriction fields), or an empty
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) if the query
returns no rows.

**Algorithm**: Filters the lookup table to the requested properties,
formats the mukey list as a quoted SQL `IN (...)` clause, builds the
`SELECT` column list (component/horizon core columns + requested
`_l/_r/_h` property columns + optional restriction columns), and joins
`mapunit` → `component` → `chorizon` → (`chfrags`, and if restrictions
requested, `corestrictions` + `chtexture`), ordered by
`cokey, hzdept_r`. Executes via
[`soilDB::SDA_query()`](http://ncss-tech.github.io/soilDB/reference/SDA_query.md);
on error, logs the failing query and re-raises via
[`stop()`](https://rdrr.io/r/base/stop.html).

### 5. **`aggregate_rock_fragment_volume_working(ssurgo_data, verbose = FALSE)`**

**Purpose**: Collapses per-fragment-size-class rock fragment volume rows
into a single summed RFV value per horizon.

**Parameters**: `ssurgo_data` - data frame containing `rfv_l/_r/_h`,
`chkey`, and `fragsize_r`; `verbose` - log aggregation stats.

**Returns**: `ssurgo_data` with one row per horizon and RFV columns
replaced by the per-`chkey` sum across fragment-size classes (rows
without `rfv_*` or `chkey`/`fragsize_r` columns are returned unchanged).

**Algorithm**: No-ops if the RFV columns or `chkey`/`fragsize_r` aren’t
present. Otherwise coerces available RFV columns to numeric, groups by
`chkey` + `fragsize_r` and takes the first value per class (removing
duplicate joins from other tables), then re-groups by `chkey` alone and
sums across fragment classes with `na.rm = TRUE`. Drops the original
(un-aggregated) RFV columns from `ssurgo_data`, left-joins the
aggregated sums back on `chkey`, and de-duplicates rows.

### 6. **`add_restriction_indicators_working(ssurgo_data, verbose = FALSE)`**

**Purpose**: Derives a suite of restriction/suitability flag columns
from horizon names, texture class, and restriction-kind fields, then
applies
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)
to produce the final `unsuitable_horizon` flag.

**Parameters**: `ssurgo_data` - data frame with restriction-related
columns (`reskind`, `resdept_l/_r`, `texcl`, `lieutex`, `hzname`,
`desgnmaster`, `hzdept_r`); `verbose` - log detection counts.

**Returns**: `ssurgo_data` with added columns
`missing_required_properties`, `restriction_top`,
`hz_below_restriction`, `has_reskind_restriction`, `is_cemented`,
`horizon_suggests_restriction`, `organic_horizon`,
`texture_suggests_restriction`, `is_potentially_restrictive`, and
`unsuitable_horizon`.

**Algorithm**: Checks which representative-value (`_r`) property columns
are actually present and flags rows missing any of them
(`missing_required_properties`); computes an effective restriction top
depth (`dplyr::coalesce(resdept_l, resdept_r)`) and whether the horizon
lies at/below it; pattern-matches `reskind` against bedrock/pan/cemented
keywords combined with the missing-properties and below-restriction
flags (`has_reskind_restriction`); flags cemented/indurated texture
(`is_cemented`); pattern-matches `hzname`/`desgnmaster` for
restriction-suggestive designations (`Cr`/`Cd`/`Cx`, “R” master,
fragipan/duripan/etc., or an `m`-suffix) and for organic horizons
(`organic_horizon`); flags cemented/indurated/bedrock texture class or
in-lieu-texture (`texture_suggests_restriction`); ORs the restriction
indicators into `is_potentially_restrictive`. Finally calls
`is_unsuitable(ssurgo_data, hzname_col = "hzname")` inside a `tryCatch`,
falling back to a simplified regex + organic/restrictive-flag
combination if that call errors.

## Utility Functions (ssurgo-acquisition.R)

- **`validate_download_inputs_ssurgo(aoi_wkt, properties, include_restrictions = TRUE, strict_geometry = TRUE, max_area_deg2 = 100)`** -
  SSURGO-specific input validator; runs
  [`validate_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_parameters.md)
  for basic type/range checks, then
  [`validate_wkt_geometry()`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
  (bounds `[-180,180] x [-90,90]`, area limits, vertex/part complexity
  limits) and
  `validate_properties_with_synonyms(property_lookup = "ssurgo")`,
  merging errors/warnings into one result with `valid`, `errors`,
  `warnings`, `corrected_parameters`, and `validation_metadata`.
- **`validate_download_inputs_ssurgo_with_config(aoi_wkt, properties, include_restrictions = TRUE, config_file = NULL, log_validation = TRUE)`** -
  config-file-aware wrapper around the same geometry/property validation
  logic, loading defaults via
  [`load_configuration()`](https://jjmaynard.github.io/soilSIM/reference/load_configuration.md)/`get_default_configuration("validation")`
  and optionally generating a report via
  [`generate_validation_report_ssurgo()`](https://jjmaynard.github.io/soilSIM/reference/generate_validation_report_ssurgo.md).
- **`generate_validation_report_ssurgo(validation_results)`** - turns a
  validation-results list into a summary report (`PASSED`/`FAILED`
  status, error/warning counts, geometry and property summaries),
  defensively handling missing/malformed metadata.
- **`check_ssurgo_cache(aoi_wkt, properties, include_restrictions, cache_dir, max_age_days = 30, verbose = FALSE)`** -
  looks up a cache file keyed by
  [`generate_ssurgo_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/generate_ssurgo_cache_key.md),
  rejecting it if missing, older than `max_age_days`, or structurally
  invalid (must contain `data`/`mu`/`metadata`).
- **`cache_ssurgo_data(data, mu, aoi_wkt, properties, include_restrictions, cache_dir, compress = TRUE, verbose = FALSE)`** -
  writes an `.rds` cache file (creating `cache_dir` if needed)
  containing the data, spatial `mu`, and metadata (timestamp, request
  parameters, row/cokey counts, R version).
- **`generate_ssurgo_cache_key(aoi_wkt, properties, include_restrictions)`** -
  builds a filename-safe cache key from truncated
  [`digest::digest()`](https://eddelbuettel.github.io/digest/man/digest.html)
  hashes of the AOI WKT and sorted property list, plus a
  restrictions-flag suffix.
- **`create_download_metadata(start_time, end_time, aoi_wkt, properties, include_restrictions, mukey_count, data_rows, unique_cokeys, spatial_result, validation_results)`** -
  assembles the timing, request, spatial-area, and validation-pass/fail
  metadata list returned alongside downloaded data.
- **`download_and_prepare_ssurgo(aoi_wkt, properties = c("clay","sand","silt","db","ph","cec","rfv","w3b","w15b"), max_depth = 250, cache_dir = NULL, verbose = FALSE)`** -
  convenience wrapper that always downloads with restrictions +
  validation on, then depth-filters `ssurgo_data` to
  `hzdepb_r <= max_depth` and attaches a `preparation_metadata` element
  (`max_depth_applied`, `ready_for_simulation = TRUE`,
  `compatible_with_infill_functions = TRUE`).

### 7. **`process_ssurgo_data()`** - Main Processing Entry Point

**Purpose**: Turns raw downloaded SSURGO data into a cleaned,
standardized, infill-ready data set, producing horizon-level,
component-level, and combined views plus a quality report.

**Parameters**:

``` r

process_ssurgo_data(
  raw_data,                    # Raw SSURGO data from download functions.
  processing_options = list(), # Overrides merged over defaults (see below).
  validate_results = TRUE,     # Logical. Run Module-8 validation on the output.
  max_depth = 250,             # Numeric. Maximum depth (cm) retained.
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

Default `processing_options` (each independently overridable):
`detect_unsuitable = TRUE`, `advanced_cleaning = TRUE`,
`standardize_names = TRUE`, `remove_invalid = TRUE`,
`calculate_derived = TRUE`, `validate_logic = TRUE`,
`preserve_compatibility = TRUE`.

**Returns**: A list with `processed_data` (the main combined output),
`horizon_data`, `component_data`, `processing_metadata` (timing, options
used, row counts, per-stage stats), `validation_results`, and
`quality_report`.

**Algorithm**: Merges user `processing_options` over the defaults via
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md),
then runs three sub-pipelines in sequence:
[`process_horizon_data_working_compatible()`](https://jjmaynard.github.io/soilSIM/reference/process_horizon_data_working_compatible.md)
for horizon-level cleaning,
[`process_component_data_working_compatible()`](https://jjmaynard.github.io/soilSIM/reference/process_component_data_working_compatible.md)
for component-level extraction, and
[`create_infill_compatible_dataset()`](https://jjmaynard.github.io/soilSIM/reference/create_infill_compatible_dataset.md)
to build the single combined data frame consumed by the infilling
functions. If `validate_results` (or `options$validate_logic`) is true,
runs
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)
and
[`generate_processing_quality_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_processing_quality_report.md)
on the result.

### 8. **`process_horizon_data_working_compatible()`**

**Purpose**: Horizon-level cleaning pipeline: standardizes names, flags
unsuitable horizons, cleans each known soil-property column, drops
invalid records, computes derived properties, and tracks statistics at
each stage.

**Parameters**:

``` r

process_horizon_data_working_compatible(
  raw_data,
  detect_unsuitable = TRUE,   # Flag unsuitable horizons via is_unsuitable().
  advanced_cleaning = TRUE,   # Run clean_property_data_ssurgo_compatible() per property.
  standardize_names = TRUE,   # Run standardize_property_names(target_standard = "ssurgo").
  remove_invalid = TRUE,      # Drop invalid horizon rows.
  calculate_derived = TRUE,   # Add hzthk_r/hz_midpoint/awc_r if derivable.
  max_depth = 250,
  verbose = FALSE
)
```

**Returns**: A list with `processed_data` (cleaned horizon data frame)
and `processing_stats` (initial/final row counts, retention rate,
unsuitable-horizon count, per-property cleaning reports, removed-invalid
count, depth range/avg thickness, unique cokey count, property
completeness).

**Algorithm**: Sequentially applies each enabled step in a fixed order
(standardize → detect-unsuitable → per-property cleaning via
[`identify_soil_property_columns_working()`](https://jjmaynard.github.io/soilSIM/reference/identify_soil_property_columns_working.md) +
[`clean_property_data_ssurgo_compatible()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data_ssurgo_compatible.md)
→
[`remove_invalid_horizons_working_compatible()`](https://jjmaynard.github.io/soilSIM/reference/remove_invalid_horizons_working_compatible.md)
→
[`calculate_derived_horizon_properties_working()`](https://jjmaynard.github.io/soilSIM/reference/calculate_derived_horizon_properties_working.md)),
recording counts/statistics into `processing_stats` after each stage,
finishing with
[`calculate_property_completeness_working()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_completeness_working.md).

### 9. **`process_component_data_working_compatible(raw_data, standardize_names = TRUE, remove_invalid = TRUE, verbose = FALSE)`**

**Purpose**: Extracts and cleans the component-level
(one-row-per-`cokey`) subset of the raw data.

**Returns**: A list with `processed_data` (distinct component rows for
the available component columns) and `processing_stats` (row counts,
retention rate, unique cokeys, component-percent range/average).

**Algorithm**: Intersects a fixed list of known component columns
(`cokey`, `compname`, `comppct_r`, taxonomic fields, `mukey`,
`compkind`, `hydgrp`, `hydric`, etc.) with what’s present in `raw_data`,
selects distinct rows on those columns, optionally standardizes names
([`standardize_property_names()`](https://jjmaynard.github.io/soilSIM/reference/standardize_property_names.md)),
removes invalid rows
([`remove_invalid_components_working()`](https://jjmaynard.github.io/soilSIM/reference/remove_invalid_components_working.md)),
and computes derived flags
([`calculate_component_stats_working()`](https://jjmaynard.github.io/soilSIM/reference/calculate_component_stats_working.md)).
Returns an empty result immediately if none of the known component
columns are present.

### 10. **`create_infill_compatible_dataset(raw_data, horizon_processing, component_processing, options, verbose = FALSE)`**

**Purpose**: Builds the single main data frame - starting from
`raw_data` rather than the separately-processed horizon/component
outputs - that downstream infilling functions consume.

**Returns**: A processed data frame: unsuitable-horizon flagged (if
`options$detect_unsuitable`), property-cleaned (if
`options$advanced_cleaning`), with essential columns guaranteed
([`ensure_essential_columns_working()`](https://jjmaynard.github.io/soilSIM/reference/ensure_essential_columns_working.md)),
depth-filtered to `options$max_depth %||% 250`, and sorted by
`cokey, hzdept_r`.

**Algorithm**: Re-derives `unsuitable_horizon` and re-runs
[`clean_property_data_ssurgo_compatible()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data_ssurgo_compatible.md)
per identified property column directly on `raw_data` (rather than
reusing `horizon_processing`/`component_processing`’s outputs, which
exist mainly for their statistics), then guarantees essential columns
exist/typed, filters `hzdepb_r <= max_depth` (keeping `NA` depths), and
arranges rows by component then top depth.

### 11. **`clean_property_data_ssurgo_compatible(df, property_name, validation_config = NULL, generate_report = TRUE, verbose = FALSE)`**

**Purpose**: Per-property cleaning: string parsing, type conversion,
infinite-value handling, statistical outlier removal, and SSURGO
range-limit enforcement, for a property’s `_l`/`_r`/`_h` columns.

**Returns**: A list with `data` (cleaned data frame) and, if
`generate_report = TRUE`, `report` (cleaning actions, type-conversion
success rates, outliers detected, string-parsing metadata).

**Algorithm**: For each available `_l/_r/_h` column: parses
character/factor values via
[`advanced_string_parser_vectorized()`](https://jjmaynard.github.io/soilSIM/reference/advanced_string_parser_vectorized.md)
or numerically converts via
[`vectorized_type_conversion()`](https://jjmaynard.github.io/soilSIM/reference/vectorized_type_conversion.md);
coerces non-finite values to `NA`; for the `_r` column only, optionally
applies caller-supplied range-rule validation
([`validate_numeric_ranges()`](https://jjmaynard.github.io/soilSIM/reference/validate_numeric_ranges.md)),
then IQR-based outlier detection via
`detect_outliers(method = "iqr", threshold = 3.0)` (outliers set to
`NA`), then
[`apply_basic_range_limits()`](https://jjmaynard.github.io/soilSIM/reference/apply_basic_range_limits.md)
for hardcoded SSURGO plausibility bounds. Note (per this function’s
`@seealso`): the string-parsing/type-conversion/range-limit helpers it
calls live in `data-infilling.R`, not this file -
`ssurgo-processing.R`’s own near-duplicate `*_working()` versions of
these helpers were consolidated into the canonical `data-infilling.R`
implementations.

### 12. **`hz_quant_prob_mukey(hz_data)`**

**Purpose**: Computes per-mukey, per-depth 5th/50th/95th percentile
statistics and 90% prediction-interval widths (`_PIW90`) across
simulated horizon data for a wide basket of soil properties, plus (when
texture data and the optional `soiltexture` package are available) the
most probable USDA texture class and its simulation-frequency
probability. This is a downstream/reporting function operating on
already-simulated (`sim_*`-style) data rather than raw acquisition
output, but lives in this file as it was ported directly from the legacy
`code/sim-functions.R`.

**Parameters**: `hz_data` - simulated horizon data frame with `mukey`,
`hzdept_r`, `hzdepb_r`, and any of a recognized property set
(sand/silt/clay totals, bulk density, `w3b`/`w15b` water retention, RFV,
pH, CEC, SOC, and optionally van Genuchten parameters
`theta_r`/`theta_s`/`alpha`/`npar`/`ksat`).

**Returns**: A data frame with one row per `mukey`/depth combination,
`_05`/`_50`/`_95`/`_PIW90`-suffixed quantile columns (alphabetically
ordered) for each available property, and, when texture columns and
`soiltexture` are present, the most probable texture class and its
probability.

**Algorithm**: Renames the available recognized property columns to
short internal names, computes grouped (`mukey`, `top`) quantiles at
`q = c(0.05, 0.5, 0.95)` via three separate
[`dplyr::summarize()`](https://dplyr.tidyverse.org/reference/summarise.html)
passes, derives `_PIW90` as the 95th-minus-5th-percentile difference,
joins all four quantile tables together and reorders columns
alphabetically, then re-attaches `bottom` depth. If `sand`/`silt`/`clay`
are present and `soiltexture` is installed, classifies each simulated
row into a USDA texture class via
[`soiltexture::TT.points.in.classes()`](https://rdrr.io/pkg/soiltexture/man/TT.points.in.classes.html),
tallies per-(mukey, top) class frequency across `simulation_number` and
keeps the modal class and its probability; otherwise emits a
[`warning()`](https://rdrr.io/r/base/warning.html) and skips texture
classification.

### 13. **`generate_processing_quality_report(original_data, processed_data, horizon_stats, component_stats, validation_results)`**

**Purpose**: Blends whichever quality signals are available (validation
score, row-retention rate, horizon property completeness, component
retention) into one weighted 0-1 score and letter grade, renormalizing
weights over only the signals actually present so a missing signal
doesn’t zero out the score.

**Returns**: A list with `overall_quality_score`, `quality_grade`
(`"Excellent"`/`"Good"`/`"Acceptable"`/`"Needs Improvement"`/`"Poor"`/`"Unknown"`),
`data_retention_rate`, a human-readable `processing_summary` string, and
the passed-through stats/validation inputs.

**Algorithm**: Coalesces each of the four candidate signals to
`NA_real_` if unavailable (validation score weight 0.40, retention 0.25,
horizon completeness 0.25, component retention 0.10), takes the weighted
mean over only the non-`NA` signals (renormalizing the weights actually
used), and maps the resulting score to a letter grade via fixed
thresholds (\>=0.9 Excellent, \>=0.8 Good, \>=0.7 Acceptable, \>=0.6
Needs Improvement, else Poor).

## Utility Functions (ssurgo-processing.R)

- **`identify_soil_property_columns_working(df)`** - returns the base
  names (e.g. `"sandtotal"`) of columns in `df` ending in `_r` that
  match a fixed list of 13 known SSURGO properties.
- **`calculate_property_completeness_working(df)`** - fraction of
  non-`NA` `_r` values per identified property.
- **`remove_invalid_horizons_working_compatible(df, max_depth = 250, verbose = FALSE)`** -
  drops rows with negative/missing/inconsistent depths
  (`hzdepb_r <= hzdept_r`, either depth `> max_depth`) or missing/empty
  `cokey`.
- **`remove_invalid_components_working(df, verbose = FALSE)`** - drops
  rows with missing/empty `cokey` or `comppct_r` outside `[0, 100]`.
- **`calculate_component_stats_working(df, verbose = FALSE)`** - adds
  `is_major_component` (`comppct_r >= 15`) and
  `has_taxonomic_classification` (non-`NA` `taxclname`) flags.
- **`calculate_derived_horizon_properties_working(df, verbose = FALSE)`** -
  adds `hzthk_r`, `hz_midpoint`, and `awc_r`
  (`wthirdbar_r - wfifteenbar_r`) when derivable and not already
  present.
- **`ensure_essential_columns_working(df)`** - synthesizes placeholder
  `cokey`/`hzname` if missing, coerces `hzdept_r`/`hzdepb_r` to numeric.

### 14. **`infill_soil_property()`** - Core Per-Property Infilling Function

**Purpose**: Fills missing representative (`_r`) and range (`_l`/`_h`)
values for a single soil property, using a six-strategy recovery
hierarchy that always excludes unsuitable horizons (R, Cr/Cd/Cx, O
horizons, etc.) from both being infilled and being used as sources.

**Parameters**:

``` r

infill_soil_property(
  df,                          # Input soil data frame.
  property_name,               # Single property base name (e.g. "claytotal"); "rfv" is special-cased.
  property_config = NULL,      # Optional config from get_default_property_config()/create_custom_property_config(); auto-looked-up if NULL.
  max_depth = 250,              # Depth (cm) constraint on which rows are eligible for infilling.
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

**Returns**: `df` with the property’s `_r` column infilled where
possible, `_l`/`_h` range values infilled, `unsuitable_horizon` and
`infill_method` (audit trail string per row) columns added/updated.

**Algorithm**: Delegates to
[`infill_rfv_property_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_rfv_property_integrated.md)
for `property_name == "rfv"`. Otherwise: ensures bookkeeping columns
exist
([`ensure_infilling_columns()`](https://jjmaynard.github.io/soilSIM/reference/ensure_infilling_columns.md)),
cleans the property
([`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md)),
computes `unsuitable_horizon` via
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md),
resolves `property_config` via
[`get_default_property_config()`](https://jjmaynard.github.io/soilSIM/reference/get_default_property_config.md)
if not supplied, and identifies “problematic” cells (missing, suitable,
within `max_depth`). If nothing is problematic, only applies range
infilling and returns early. Otherwise groups by `cokey` (or `compname`,
via
[`determine_grouping_column()`](https://jjmaynard.github.io/soilSIM/reference/determine_grouping_column.md))
and applies, in order: **Strategy 1-3** inside
[`process_property_group()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group.md)
→
[`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md)
(horizon-name matching, then depth-weighted averaging, then
within-component interpolation - each only proceeding if the prior
strategy left cells unfilled); **Strategy 4**
[`cross_component_property_interpolation()`](https://jjmaynard.github.io/soilSIM/reference/cross_component_property_interpolation.md)
(whole-dataset level, uses other components’ data at similar depths);
**Strategy 5**
[`related_property_estimation()`](https://jjmaynard.github.io/soilSIM/reference/related_property_estimation.md)
(whole-dataset level, pedological relationships e.g. texture-sum,
clay/OM-based CEC); **Strategy 6**
[`process_property_group_fallback()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group_fallback.md)
→
[`apply_group_fallback_mean()`](https://jjmaynard.github.io/soilSIM/reference/apply_group_fallback_mean.md)
(per-group depth-weighted or plain mean, last resort). Finishes with
[`infill_property_range_values()`](https://jjmaynard.github.io/soilSIM/reference/infill_property_range_values.md)
to fill `_l`/`_h` and enforce `_l <= _r <= _h` ordering.

### 15. **`process_soil_properties_comprehensive()`** - Multi-Property Workflow

**Purpose**: Orchestrates
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)
(and RFV’s special path) across a whole property set in three phases,
with optional post-hoc filtering.

**Parameters**:

``` r

process_soil_properties_comprehensive(
  df,
  properties = NULL,           # NULL = auto-detect via auto_detect_soil_properties().
  max_depth = 250,
  remove_unsuitable = FALSE,   # Drop unsuitable-horizon rows from the output.
  remove_incomplete = FALSE,   # Drop rows incomplete in required_properties.
  required_properties = NULL,  # Defaults to `properties` if remove_incomplete and NULL.
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

**Returns**: `df` with all requested properties infilled, plus
`processing_log`, `processing_time`, and `properties_processed`
attributes attached (`attr(df, ...)`).

**Algorithm**: Auto-detects properties if not given
([`auto_detect_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/auto_detect_soil_properties.md)),
validates them
([`validate_properties_with_synonyms()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties_with_synonyms.md),
warnings only), ensures bookkeeping columns
([`ensure_infilling_columns()`](https://jjmaynard.github.io/soilSIM/reference/ensure_infilling_columns.md)),
then processes in three phases: **Phase 1** foundation properties
(`sandtotal`, `claytotal`, `silttotal`, `dbovendry`, `rfv`) via
[`process_single_property()`](https://jjmaynard.github.io/soilSIM/reference/process_single_property.md);
**Phase 2** water retention (`wthirdbar`, `wfifteenbar`) via
[`infill_water_retention_saxton_rawls_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_water_retention_saxton_rawls_integrated.md),
only if texture+bulk-density columns are all present; **Phase 3**
remaining (“chemical”) properties via
[`process_single_property()`](https://jjmaynard.github.io/soilSIM/reference/process_single_property.md).
Optionally applies
[`apply_data_filtering()`](https://jjmaynard.github.io/soilSIM/reference/apply_data_filtering.md)
afterward, and logs a summary
([`generate_processing_summary()`](https://jjmaynard.github.io/soilSIM/reference/generate_processing_summary.md))
if verbose.

### 16. **`infill_rfv_property_integrated(df, max_depth = 250, verbose = FALSE)`**

**Purpose**: Specialized rock-fragment-volume infilling, applying
row-level texture-informed imputation rather than the general
multi-strategy hierarchy.

**Returns**: `df` with `rfv_l/_r/_h` populated for suitable, in-depth
rows, and `infill_method` annotated `"rfv_specialized; "`.

**Algorithm**: Determines the suitable + within-`max_depth` row set,
then for each such row calls
[`impute_rfv_values()`](https://jjmaynard.github.io/soilSIM/reference/impute_rfv_values.md)
(row-level: texture-sum-based default of 0.5% if texture sums to ~100
and RFV missing, else a generic 2% default; near-zero values floored to
0.1%; valid values get a +/-30% `_l/_h` spread clamped to `[0.05, 85]`).
Finishes with
[`infill_property_range_values()`](https://jjmaynard.github.io/soilSIM/reference/infill_property_range_values.md)
for any remaining range gaps.

### 17. **`infill_water_retention_saxton_rawls_integrated()`**

**Purpose**: Estimates missing field-capacity (`wthirdbar`) and
wilting-point (`wfifteenbar`) values from texture, bulk density, rock
fragments, and organic matter via the Saxton-Rawls pedotransfer
function.

**Parameters**:

``` r

infill_water_retention_saxton_rawls_integrated(
  df,
  max_depth = 250,
  add_ranges = TRUE,    # Also populate _l/_h range columns.
  overwrite = FALSE,    # If TRUE, re-estimate even where _r already has a value.
  verbose = FALSE
)
```

**Returns**: `df` with `wthirdbar_r`/`wfifteenbar_r` (and `_l`/`_h` if
`add_ranges`) filled where texture + bulk density inputs are complete;
unchanged if `claytotal_r`/`sandtotal_r`/`silttotal_r`/`dbovendry_r`
aren’t all present.

**Algorithm**: Requires all four texture/bulk-density columns to be
present (returns `df` unchanged with a warning log otherwise).
Computes/reuses `unsuitable_horizon`, builds a suitable+in-depth
processing mask, and (unless `overwrite`) restricts to rows where the
target `_r` value is still missing. Row by row, calls
[`calculate_saxton_rawls_single()`](https://jjmaynard.github.io/soilSIM/reference/calculate_saxton_rawls_single.md)
with the row’s texture/bulk-density/RFV(default 0)/OM(default 2%) values
and writes back `field_capacity`/`wilting_point` (and their `_l`/`_h`
spreads if `add_ranges`), annotating `infill_method`.

## Utility Functions (data-infilling.R)

- **`clean_property_data(df, property_name, validation_config = NULL, generate_report = FALSE, verbose = FALSE)`** -
  the
  [`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)-side
  analog of
  [`clean_property_data_ssurgo_compatible()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data_ssurgo_compatible.md):
  string parsing/type conversion/infinite handling, plus soil-aware
  outlier detection
  ([`detect_statistical_outliers_soil_aware()`](https://jjmaynard.github.io/soilSIM/reference/detect_statistical_outliers_soil_aware.md))
  and
  [`apply_basic_range_limits()`](https://jjmaynard.github.io/soilSIM/reference/apply_basic_range_limits.md)
  on the `_r` column.
- **`get_default_property_config(property_name)`** - returns a built-in
  config (`type`, `units`, `typical_range`, `fallback_range`, and
  property-specific flags like
  `clay_dependent`/`horizon_effects`/`related_properties`) for texture,
  bulk density, water retention, RFV, CEC, pH, and organic-matter/carbon
  properties, or a generic fallback.
- **`create_custom_property_config(property_name, property_type = "generic", units = "unknown", typical_range = NULL, fallback_range = 5, related_properties = NULL, special_options = NULL)`** -
  builds a custom property config with input validation.
- **`validate_property_config(config, property_name)`** - checks a
  config has `type`/`units`/`fallback_range` and well-formed
  `typical_range`; stops with an error otherwise.
- **`get_rfv_range_category(rfv_value)`** - categorizes an RFV
  percentage into
  `"none"`/`"low"`/`"moderate"`/`"high"`/`"very_high"`/`"extreme"`.
- **`apply_property_constraints(values, property_config)`** - clamps
  values to `typical_range` and type-specific bounds (e.g. `[0,100]`
  texture, `[0,14]` pH, `[0.01,95]` rock fragments).
- **[`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md)
  /
  [`add_range_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_range_rule.md)
  /
  [`add_relationship_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_relationship_rule.md)
  /
  [`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md)** -
  a small pluggable rule-configuration system for range/relationship
  checks; note relationship rules are recorded but **not** currently
  enforced by
  [`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md)
  (matches legacy behavior, not a new bug).
- **`summarize_unsuitable_horizons(df, hzname_col = "hzname")`** -
  reporting companion to
  [`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md);
  returns count and unique horizon names of excluded horizons, logging
  via
  [`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md).
- **`infill_property_range_values(df, property_name, property_config)`** -
  fills missing `_l`/`_h` from learned
  ([`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md))
  and contextual
  ([`get_property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/get_property_contextual_ranges.md))
  spreads via
  [`calculate_property_lower_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_lower_bound.md)/[`calculate_property_upper_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_upper_bound.md),
  then enforces `_l <= _r <= _h`.
- **`learn_property_ranges(df, property_name, property_config)`** -
  learns median `_r - _l` / `_h - _r` spreads from complete,
  suitable-horizon rows, broken out by horizon name, by depth zone
  (surface/subsurface/deep), and overall.
- **`get_property_contextual_ranges(df, property_name, property_config)`** -
  hardcoded pedological-knowledge spread tables per property type
  (texture, bulk density by horizon letter, water retention, RFV
  category).
- **[`calculate_property_lower_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_lower_bound.md)
  /
  [`calculate_property_upper_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_upper_bound.md)** -
  per-row bound calculation using a priority order: horizon-learned →
  depth-learned → overall-learned → contextual → fallback spread.
- **`get_contextual_spread(row, property_name, context_ranges, property_config, bound_type)`** -
  looks up the appropriate contextual spread for
  [`calculate_property_lower_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_lower_bound.md)/`_upper_bound()`.
- **`auto_detect_soil_properties(df)`** - returns which of a fixed
  physical/chemical/water-retention property list have a corresponding
  `_r` column present.
- **`ensure_infilling_columns(df, verbose = FALSE)`** - guarantees
  `infill_method`, `cokey` (synthesized from `compname` or row index if
  absent), and `hzname` (`"Unknown"` if absent) exist; coerces depth
  columns numeric.
- **`determine_grouping_column(df)`** - returns `"cokey"` if present,
  else `"compname"`, else `NULL`.
- **[`process_single_property()`](https://jjmaynard.github.io/soilSIM/reference/process_single_property.md)
  /
  [`apply_data_filtering()`](https://jjmaynard.github.io/soilSIM/reference/apply_data_filtering.md)
  /
  [`generate_processing_summary()`](https://jjmaynard.github.io/soilSIM/reference/generate_processing_summary.md)** -
  per-property logging wrapper, post-hoc unsuitable/incomplete-row
  filtering, and a verbose text summary, all used by
  [`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md).
- **`advanced_string_parser_vectorized(value_strings)` /
  `parse_single_string_advanced(value_string)`** - the canonical string
  parser (ranges like `"15-20"`, comparison operators `"<5"`/`">50"`, a
  qualitative-term table like `"TRACE"`/`"LOW"`/`"HIGH"`, percentages,
  and plain numerics), each parse tagged with a `type` and `confidence`
  score.
- **`vectorized_type_conversion(values)`** -
  factor/character/numeric-safe coercion to numeric, stripping
  non-numeric characters from strings.
- **`detect_statistical_outliers_soil_aware(values, df, property_name)`** -
  property-aware outlier rules: conservative absolute bounds for
  pH/texture/bulk-density (only physically-impossible values flagged),
  else a very-conservative (factor 5.0) IQR rule.
- **`apply_basic_range_limits(values, property_name)`** - hardcoded
  per-property plausibility bounds (texture `[0,100]`, bulk density
  `[0.3,3.0]`, water retention, CEC `[0,200]`, pH
  `[2.5,11.0]`/`[2.0,10.5]`, OM/OC, RFV `[0,95]`); values outside are
  set `NA`; unknown properties get a non-negative floor only.
- **`standardize_horizon_name(hzname)` /
  `calculate_horizon_similarity(hz1, hz2)`** - name normalization
  (uppercase, strip trailing digits/punctuation) and a 0-1 similarity
  score (exact match = 1.0, same leading letter = 0.8 +
  character-overlap bonus, related horizon-letter groups
  e.g. `A`/`AP`/`AE` = 0.6, else 0).
- **`impute_rfv_values(row)`** - the row-level RFV imputation logic used
  by
  [`infill_rfv_property_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_rfv_property_integrated.md).
- **`calculate_saxton_rawls_single(sand_pct, clay_pct, silt_pct, bulk_density, rfv_pct = 0, om_pct = 2)`** -
  the Saxton-Rawls pedotransfer math itself (clamps inputs, renormalizes
  texture to 100 if off by \>5 points, computes saturated water
  content/field-capacity/wilting-point via the published regression
  equations, applies bulk-density and RFV volumetric correction, and
  returns `field_capacity`/`wilting_point` plus `_l`/`_h` (+/-15%) and
  `available_water_capacity`).
- **[`depth_weighted_property_infill()`](https://jjmaynard.github.io/soilSIM/reference/depth_weighted_property_infill.md)
  /
  [`within_component_property_interpolation()`](https://jjmaynard.github.io/soilSIM/reference/within_component_property_interpolation.md)** -
  Strategies 2 and 3 inside
  [`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md):
  nearest-depth (\<=20cm tolerance) inverse-distance-weighted mean, then
  [`approx()`](https://rdrr.io/r/stats/approxfun.html)-based linear
  interpolation within a component.
- **`cross_component_property_interpolation(group, property_col)`** -
  Strategy 4: fills missing suitable-horizon values using other
  components’ suitable horizons within 15cm depth tolerance,
  inverse-distance weighted.
- **`related_property_estimation(group, property_name, property_config)`** -
  Strategy 5: property-type-specific pedological relationships (texture
  sum-to-100; clay-based water retention; clay+OM-based CEC;
  horizon/OM-adjusted pH; depth/horizon/clay-adjusted OM;
  texture-adjusted bulk density).
- **`calculate_depth_weighted_mean(group, property_col)`** /
  **[`apply_group_fallback_mean()`](https://jjmaynard.github.io/soilSIM/reference/apply_group_fallback_mean.md)** -
  Strategy 6 (last resort): thickness-weighted (or plain) mean of
  suitable horizons in a group.
- **[`horizon_name_property_infill()`](https://jjmaynard.github.io/soilSIM/reference/horizon_name_property_infill.md)
  /
  [`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md)** -
  Strategy 1 (horizon-name similarity matching, weighted mean of matches
  with similarity \> 0.5) and its orchestration across Strategies 1-3.
- **[`process_property_group()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group.md)
  /
  [`process_property_group_fallback()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group_fallback.md)** -
  per-group
  ([`dplyr::group_modify()`](https://dplyr.tidyverse.org/reference/group_map.html)-invoked)
  wrappers that recompute the problematic-cell mask locally before
  delegating to
  [`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md)
  /
  [`apply_group_fallback_mean()`](https://jjmaynard.github.io/soilSIM/reference/apply_group_fallback_mean.md)
  respectively.

## Internal Connections

    download_ssurgo_tabular() [ssurgo-acquisition.R, MASTER]
    ├── validate_download_inputs_ssurgo()
    ├── check_ssurgo_cache()                       [cache hit -> early return]
    ├── create_ssurgo_property_lookup_working()
    ├── process_aoi_and_get_mukeys_working()
    ├── execute_ssurgo_query_working()
    ├── aggregate_rock_fragment_volume_working()    [if "rfv" requested]
    ├── add_restriction_indicators_working()        [if include_restrictions]
    │   └── is_unsuitable()                         [utils.R]
    ├── validate_data_quality()                     [utils.R]
    ├── cache_ssurgo_data()
    │   └── generate_ssurgo_cache_key()
    └── create_download_metadata()

    download_and_prepare_ssurgo() [ssurgo-acquisition.R, convenience wrapper]
    └── download_ssurgo_tabular()

    process_ssurgo_data() [ssurgo-processing.R, MAIN ENTRY]
    ├── process_horizon_data_working_compatible()
    │   ├── standardize_property_names()            [utils.R]
    │   ├── is_unsuitable()                         [utils.R]
    │   ├── identify_soil_property_columns_working()
    │   ├── clean_property_data_ssurgo_compatible()  [per property]
    │   │   ├── advanced_string_parser_vectorized()  [data-infilling.R]
    │   │   ├── vectorized_type_conversion()         [data-infilling.R]
    │   │   ├── detect_outliers()                    [utils.R]
    │   │   └── apply_basic_range_limits()           [data-infilling.R]
    │   ├── remove_invalid_horizons_working_compatible()
    │   ├── calculate_derived_horizon_properties_working()
    │   └── calculate_property_completeness_working()
    ├── process_component_data_working_compatible()
    │   ├── standardize_property_names()             [utils.R]
    │   ├── remove_invalid_components_working()
    │   └── calculate_component_stats_working()
    ├── create_infill_compatible_dataset()
    │   ├── is_unsuitable()                          [utils.R]
    │   ├── clean_property_data_ssurgo_compatible()  [per property]
    │   └── ensure_essential_columns_working()
    ├── validate_data_quality()                      [utils.R]
    └── generate_processing_quality_report()

    hz_quant_prob_mukey() [ssurgo-processing.R, standalone - downstream/reporting use]

    infill_soil_property() [data-infilling.R, PER-PROPERTY CORE]
    ├── infill_rfv_property_integrated()             [if property_name == "rfv"]
    │   └── impute_rfv_values()
    ├── ensure_infilling_columns()
    ├── clean_property_data()
    │   ├── advanced_string_parser_vectorized() / vectorized_type_conversion()
    │   ├── detect_statistical_outliers_soil_aware()
    │   └── apply_basic_range_limits()
    ├── is_unsuitable()                              [utils.R]
    ├── get_default_property_config()                [if property_config = NULL]
    ├── determine_grouping_column()
    ├── process_property_group()                     [Strategies 1-3, per group]
    │   └── infill_missing_property_data()
    │       ├── horizon_name_property_infill()       [Strategy 1]
    │       ├── depth_weighted_property_infill()      [Strategy 2]
    │       └── within_component_property_interpolation() [Strategy 3]
    ├── cross_component_property_interpolation()      [Strategy 4, whole dataset]
    ├── related_property_estimation()                 [Strategy 5, whole dataset]
    ├── process_property_group_fallback()             [Strategy 6, per group]
    │   └── apply_group_fallback_mean()
    │       └── calculate_depth_weighted_mean()
    └── infill_property_range_values()
        ├── learn_property_ranges()
        ├── get_property_contextual_ranges()
        ├── calculate_property_lower_bound()
        │   └── get_contextual_spread()
        └── calculate_property_upper_bound()
            └── get_contextual_spread()

    process_soil_properties_comprehensive() [data-infilling.R, MULTI-PROPERTY WORKFLOW]
    ├── auto_detect_soil_properties()                [if properties = NULL]
    ├── validate_properties_with_synonyms()          [utils.R]
    ├── ensure_infilling_columns()
    ├── process_single_property()                    [Phase 1 + Phase 3, per property]
    │   ├── infill_rfv_property_integrated()         [property == "rfv"]
    │   └── infill_soil_property()                   [all other properties]
    ├── infill_water_retention_saxton_rawls_integrated() [Phase 2]
    │   └── calculate_saxton_rawls_single()          [per row]
    ├── apply_data_filtering()                       [if remove_unsuitable/remove_incomplete]
    └── generate_processing_summary()                [if verbose]

## Dependencies

**From elsewhere in soilSIM** (the “Module 8”-equivalent utility layer,
all in `R/utils.R`):
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md),
[`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md)
*(referenced in this group’s older code comments; not called directly in
the current source)*,
[`validate_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_parameters.md),
[`validate_wkt_geometry()`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md),
[`validate_properties_with_synonyms()`](https://jjmaynard.github.io/soilSIM/reference/validate_properties_with_synonyms.md),
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
[`validate_numeric_ranges()`](https://jjmaynard.github.io/soilSIM/reference/validate_numeric_ranges.md),
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md),
[`standardize_property_names()`](https://jjmaynard.github.io/soilSIM/reference/standardize_property_names.md),
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md),
[`load_configuration()`](https://jjmaynard.github.io/soilSIM/reference/load_configuration.md),
[`get_default_configuration()`](https://jjmaynard.github.io/soilSIM/reference/get_default_configuration.md),
[`setup_logging()`](https://jjmaynard.github.io/soilSIM/reference/setup_logging.md),
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md),
and the `%||%` null-coalescing helper.

**External packages**: `soilDB` (`mukey.wcs()` for spatial mukey lookup,
`SDA_query()` for the SDA SQL query itself); `terra` (`vect()`,
`project()`, `as.polygons()`, `ext()`, `values()`, `expanse()` for AOI
geometry handling); `dplyr` (`mutate`, `filter`, `group_by`/`summarise`,
`case_when`, joins, `across` - used throughout for data-frame
transformations); `stringr` (`str_replace`/`str_replace_all` in
[`standardize_horizon_name()`](https://jjmaynard.github.io/soilSIM/reference/standardize_horizon_name.md));
`rlang` (`sym()` for programmatic grouping); `digest` (cache key
hashing); `stats` (`quantile`, `approx`, `weighted.mean`); `soiltexture`
(optional, `Suggests`; USDA texture classification in
[`hz_quant_prob_mukey()`](https://jjmaynard.github.io/soilSIM/reference/hz_quant_prob_mukey.md),
guarded by [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)).

**Downstream consumers**: The cleaned/infilled horizon+component data
frames (with complete `_l/_r/_h` triplet columns) produced by this group
feed the **Statistics & Diagnostics** module (descriptive statistics,
correlation analysis over the property triplets) and the **Monte Carlo
simulation core** (percentile-triplet distribution fitting and
correlated sampling consume the `sim_*`/triplet columns this group
guarantees are gap-free).
[`hz_quant_prob_mukey()`](https://jjmaynard.github.io/soilSIM/reference/hz_quant_prob_mukey.md)
runs in the opposite direction - it consumes *already-simulated* horizon
data (post Monte Carlo) to produce per-mukey/depth summary statistics,
so it is better understood as a reporting utility co-located here for
historical (ported-file) reasons rather than a pure
acquisition/processing step.

## Data Flow In/Out

**Inputs required from the caller**: - An **AOI** as a WKT string in
`EPSG:4326` (geographic WGS84) coordinates, valid per
[`validate_wkt_geometry()`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(bounded to `[-180,180] x [-90,90]`, default max area 100 deg^2, max
1000 vertices / 100 parts). - A **property name list** drawn from (or
resolvable via synonym to) SSURGO base names: `sandtotal`, `claytotal`,
`silttotal`, `dbovendry`, `ph1to1h2o`, `cec7`, `om`, `wthirdbar`,
`wfifteenbar` (default set), plus `rfv` for rock fragments, which is
handled specially throughout. - Optionally, a `cache_dir` for download
caching, and a `max_depth` (cm) governing how deep into the profile
processing/infilling extends (default 250 cm throughout).

**Outputs produced**: -
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
→ raw `ssurgo_data` (one row per horizon per component, `_l/_r/_h`
columns per requested property, restriction/unsuitable-horizon flags if
requested) + spatial `mu` + metadata. -
[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)
→ `processed_data` (the infill-ready combined data frame - cleaned,
standardized column names, essential columns guaranteed, depth-filtered,
sorted by `cokey`/depth), plus separate `horizon_data`/`component_data`
views and a quality report. -
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)
/
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md)
→ the same data frame shape with **no remaining `NA`s** in suitable,
in-depth `_r` cells for the processed properties (to the extent
recoverable), complete `_l`/`_h` ranges satisfying `_l <= _r <= _h`, an
`unsuitable_horizon` logical flag, and an `infill_method` audit-trail
string per row documenting which strategy filled each cell.

A caller wiring this group into a larger pipeline should expect to
supply only the WKT + property list, and should treat
`unsuitable_horizon == TRUE` rows as excluded from simulation inputs by
convention (they are never used as infilling *sources*, and callers may
choose to filter them out entirely via `remove_unsuitable = TRUE` in
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md)).

## Usage Example

``` r

# 1. Download raw SSURGO tabular data for an AOI
dl <- download_ssurgo_tabular(
  aoi_wkt = "POLYGON((-120.5 46.5, -120.4 46.5, -120.4 46.6, -120.5 46.6, -120.5 46.5))",
  properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "om", "wthirdbar", "wfifteenbar"),
  include_restrictions = TRUE,
  cache_dir = "cache/ssurgo",
  verbose = TRUE
)

# 2. Clean and standardize into an infill-ready data set
proc <- process_ssurgo_data(
  raw_data = dl$ssurgo_data,
  max_depth = 200,
  verbose = TRUE
)

# 3. Infill missing values property by property (or all at once)
infilled <- process_soil_properties_comprehensive(
  df = proc$processed_data,
  properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "om", "wthirdbar", "wfifteenbar"),
  max_depth = 200,
  remove_unsuitable = FALSE,
  verbose = TRUE
)

# `infilled` now has complete sandtotal_r/_l/_h, claytotal_r/_l/_h, ... triplets
# ready for Statistics & Diagnostics and the Monte Carlo simulation core.
```
