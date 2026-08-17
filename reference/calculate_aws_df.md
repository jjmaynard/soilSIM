# Calculate Available Water Storage by Depth Interval

Runs
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)
on soil texture/bulk-density/water-retention inputs to derive van
Genuchten pedotransfer parameters, Monte Carlo-simulates available water
holding capacity per horizon via
[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md),
then depth-slices and summarizes mean AWHC per component over standard
depth intervals via
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html).

## Usage

``` r
calculate_aws_df(sim_data_df)
```

## Arguments

- sim_data_df:

  A data frame with one row per horizon, with columns `sand_total`,
  `silt_total`, `clay_total`, `bulk_density_third_bar`,
  `water_retention_third_bar`, `water_retention_15_bar` (ROSETTA's
  expected input variable names - note these differ from the rest of
  `soilSIM`'s `sandtotal_r`/`claytotal_r`/etc. SSURGO-derived naming
  convention, since they're ROSETTA's own API contract), plus
  `compname`, `hzdept_r`, `hzdepb_r`, `cokey`.

## Value

A long-format data frame with one row per `cokey` per depth slab
actually spanned by that component's horizons (slab boundaries per
`slab.structure = c(0, 5, 15, 30, 60, 100)`), with columns `cokey`,
`top`, `bottom`, `AWHC` (mean available water holding capacity over that
slab). The
[`tidyr::pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
step widens on `variable` (always `"AWHC"` here, since only one property
is slabbed), so it does not collapse depth slabs into columns - it
exists to match
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)'s
long output shape to a plain `value`-column contract.

## Known limitation

This function requires **live network access** -
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)
POSTs to `https://www.handbook60.org/api/v1/rosetta/<version>` (via
[`httr::POST()`](https://httr.r-lib.org/reference/POST.html)) rather
than computing pedotransfer parameters locally. Calling this without
internet connectivity, or if `handbook60.org` is unreachable, will fail.
