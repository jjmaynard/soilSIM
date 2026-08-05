# Download and Prepare SSURGO Data (Workflow Convenience Wrapper)

Thin convenience wrapper around
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
that always requests restriction data and validation, then filters the
result to `max_depth` - a one-call entry point for the common "download,
ready for simulation" path.

## Usage

``` r
download_and_prepare_ssurgo(
  aoi_wkt,
  properties = c("clay", "sand", "silt", "db", "ph", "cec", "rfv", "w3b", "w15b"),
  max_depth = 250,
  cache_dir = NULL,
  verbose = FALSE
)
```

## Arguments

- aoi_wkt:

  Area of interest in WKT format.

- properties:

  Character vector of SSURGO properties to download.

- max_depth:

  Maximum depth (cm); horizons deeper than this are dropped. Pass `NULL`
  to skip depth filtering.

- cache_dir:

  Optional cache directory (see
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)).

- verbose:

  Logical; provide progress messages.

## Value

The same list
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
returns, with `ssurgo_data` depth-filtered and a `preparation_metadata`
element added.
