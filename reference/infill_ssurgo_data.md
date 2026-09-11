# Infill Missing Soil Property Values (raster-fusion entry point)

Thin wrapper: infills the standard SSURGO property set present in `df`
via
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md) -
the single infilling orchestrator - restricted to the columns actually
present. Kept as a named entry point because the raster-fusion pipeline
([`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md))
calls it with a fixed property vocabulary.

## Usage

``` r
infill_ssurgo_data(df, water_retention_method = c("saxton_rawls", "generic"))
```

## Arguments

- df:

  A horizon data frame, as returned by
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md).

- water_retention_method:

  Passed through to
  [`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md):
  `"saxton_rawls"` (default) derives `wthirdbar`/`wfifteenbar` via the
  pedotransfer function where texture + bulk density allow; `"generic"`
  sends them through the six-strategy hierarchy.

## Value

`df` with missing values infilled where possible.
