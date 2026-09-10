# Attach OSD-Derived Boundary Distinctness to Horizon Data (per-genhz `bound_sd`)

[`query_osd_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/query_osd_distinctness.md)'s
`bound_sd` (OSD boundary-distinctness offset: Abrupt/Clear/
Gradual/Diffuse converted to a numeric scale via
[`aqp::hzDistinctnessCodeToOffset()`](https://ncss-tech.github.io/aqp/reference/hzDistinctnessCodeToOffset.html))
is also used by
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)
as
[`aqp::perturb()`](https://ncss-tech.github.io/aqp/reference/perturb.html)'s
`boundary.attr` for horizon-depth perturbation. This function exposes
the same genhz-averaged `bound_sd` value on arbitrary horizon data (via
the same
`dplyr::group_by(id, genhz) |> dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE))`
pattern), so it can also gate the vertical-correlation depth kernel.

## Usage

``` r
attach_osd_boundary_distinctness(hz_data)
```

## Arguments

- hz_data:

  A data frame with at least `compname`, `hzname`, and `genhz` columns
  (the same shape
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  builds right after
  [`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)).

## Value

`hz_data` with a `bound_sd` numeric column added/overwritten (one value
per `compname`/`genhz` combination, `NA` for any combination the OSD
lookup couldn't resolve). On OSD lookup failure (e.g. no network
access -
[`query_osd_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/query_osd_distinctness.md)
calls
[`soilDB::fetchOSD()`](http://ncss-tech.github.io/soilDB/reference/fetchOSD.md)),
degrades to an all-`NA` `bound_sd` column rather than erroring, so
callers (and the kernel gating, which treats missing distinctness data
as "no gating - fall back to the plain kernel") keep working without OSD
access.
