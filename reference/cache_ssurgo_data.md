# Cache SSURGO Data (Enhanced)

Enhanced data caching with comprehensive metadata.

## Usage

``` r
cache_ssurgo_data(
  data,
  mu,
  aoi_wkt,
  properties,
  include_restrictions,
  cache_dir,
  compress = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  Downloaded SSURGO data

- mu:

  Spatial map unit data

- aoi_wkt:

  Original AOI WKT

- properties:

  Requested properties

- include_restrictions:

  Whether restrictions were included

- cache_dir:

  Cache directory path

- compress:

  Logical; compress cached data

- verbose:

  Logical; provide progress messages

## Value

List with cache operation results
