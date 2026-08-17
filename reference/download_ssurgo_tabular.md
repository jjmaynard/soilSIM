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
