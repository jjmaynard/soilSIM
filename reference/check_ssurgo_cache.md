# Check SSURGO Data Cache (Enhanced)

Enhanced cache checking with metadata validation.

## Usage

``` r
check_ssurgo_cache(
  aoi_wkt,
  properties,
  include_restrictions,
  cache_dir,
  max_age_days = 30,
  verbose = FALSE
)
```

## Arguments

- aoi_wkt:

  Area of interest WKT

- properties:

  Vector of properties

- include_restrictions:

  Logical for restriction data

- cache_dir:

  Cache directory path

- max_age_days:

  Maximum age of cached data in days

- verbose:

  Logical; provide progress messages

## Value

List with cache check results or NULL if no valid cache
