# Query SSURGO Database for Soil Data by Mukey

Queries the NRCS SSURGO database
(`mapunit`/`component`/`chorizon`/`chfrags`) for horizon morphology,
texture, bulk density, water-retention, and rock-fragment-volume data
for the given map unit key(s).

## Usage

``` r
get_aws_data_by_mukey(mukeys)
```

## Arguments

- mukeys:

  A vector of string or numeric map unit keys (mukeys) representing the
  map units to retrieve data for.

## Value

A data frame with soil horizon data for the given mukeys, including
aggregated rock fragment data.

## Details

This deliberately does not reuse
[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)
(which also joins `mapunit`/`component`/`chorizon`/`chfrags`): that
function only ever selects `chf.fragsize_r`, never `chf.fragvol_l/r/h`,
so its rock fragment aggregation path
([`aggregate_rock_fragment_volume_working()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_rock_fragment_volume_working.md))
currently produces no data. This function's entire purpose is summing
`fragvol_l/r/h` per horizon, so delegating would silently drop rock
fragment volume - a real regression, not a style difference.
