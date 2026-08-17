# Download SSURGO Tabular Data with Comprehensive Processing

Enhanced version of the proven download_ssurgo_tabular function with
additional caching, validation, and reporting capabilities while
maintaining compatibility with existing soil property simulation
workflows.

## Usage

``` r
download_ssurgo_tabular(
  aoi_wkt,
  properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "ph1to1h2o", "cec7",
    "om", "wthirdbar", "wfifteenbar"),
  include_restrictions = TRUE,
  cache_dir = NULL,
  force_download = FALSE,
  validate_data = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- aoi_wkt:

  Character. Well-Known Text (WKT) representation of the area of
  interest

- properties:

  Character vector. Soil properties to download. Default: c("w3b",
  "w15b", "db", "cec", "rfv", "clay", "ph", "sand", "silt", "soc")

- include_restrictions:

  Logical. Whether to include horizon restriction data for unsuitable
  horizon detection (default: TRUE)

- cache_dir:

  Character. Directory for caching downloaded data (default: NULL for no
  caching)

- force_download:

  Logical. If TRUE, bypass cache and download fresh data (default:
  FALSE)

- validate_data:

  Logical. Perform comprehensive data validation (default: TRUE)

- verbose:

  Logical. Provide detailed progress messages (default: FALSE)

## Value

List containing:

- ssurgo_data: Combined horizon and component data with restriction
  flags

- mu: Spatial map unit data (from soilDB::mukey.wcs)

- metadata: Download and processing metadata

- validation_results: Data validation results (if validate_data = TRUE)

- cache_info: Cache usage information (if caching enabled)

- components_missing_horizons: Data frame
  (`mukey`/`cokey`/`compname`/`comppct_l/r/h`), one row per component
  that has a real `comppct` but zero `chorizon` rows in SDA AND no AOI
  sibling (same `compname`) to recover a profile from - empty if none.
  Always empty for a cache hit (see Component recovery section below).

- components_recovered: Data frame
  (`mukey`/`cokey`/`compname`/`n_sibling_cokeys_used`), one row per
  component whose horizon profile WAS synthesized from AOI siblings -
  empty if none. Always empty for a cache hit (see Component recovery
  section below).

## Component recovery (on by default)

[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)'s
`INNER JOIN chorizon` silently drops any real component with zero
`chorizon` rows in SDA - a common, verified SSURGO data-completeness
gap, especially for minor components. This function recovers such
components automatically: for each one, it searches this same AOI's own
fetched data for other cokeys sharing the missing component's `compname`
that DO have full horizon data, and averages their profiles (aligned by
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
group) to synthesize a representative one - see
[`synthesize_component_horizons_from_siblings()`](https://jjmaynard.github.io/soilSIM/reference/synthesize_component_horizons_from_siblings.md)'s
docs for the exact averaging rule. Synthesized rows are tagged
`component_synthesized = TRUE` and
`infill_method = "component_aoi_average"` in `ssurgo_data` for
traceability. A component with no AOI sibling to recover from is left
out of `ssurgo_data` (same as before this feature existed) but reported
in `components_missing_horizons` instead of vanishing with zero trace.

## Cache invalidation

A cache entry (`cache_dir`) written before this feature shipped, or
before it was extended, predates component recovery entirely -
`ssurgo_data` in that cached entry won't contain any synthesized rows,
and a cache hit always returns empty `components_missing_horizons`/
`components_recovered` placeholders regardless of what a fresh download
would find (only counts, not the full data frames, are persisted via
`metadata`). Clear the cache to pick up recovered components.
