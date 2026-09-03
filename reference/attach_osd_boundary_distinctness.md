# Attach OSD-Derived Boundary Distinctness to Horizon Data (per-genhz `bound_sd`)

[`query_osd_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/query_osd_distinctness.md)'s
`bound_sd` (OSD boundary-distinctness offset, e.g. Abrupt/Clear/
Gradual/Diffuse converted to a numeric scale via
[`aqp::hzDistinctnessCodeToOffset()`](https://ncss-tech.github.io/aqp/reference/hzDistinctnessCodeToOffset.html))
was previously computed and consumed ONLY inside
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md),
as
[`aqp::perturb()`](https://ncss-tech.github.io/aqp/reference/perturb.html)'s
`boundary.attr` for horizon-DEPTH perturbation - it never reached the
property-value Monte Carlo simulation pipeline
([`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md),
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)).
This exposes the same already-computed, already-genhz-averaged
`bound_sd` value on arbitrary horizon data (via the identical
`dplyr::group_by(id, genhz) |> dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE))`
pattern
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)
already uses), so it can also gate the vertical-correlation depth kernel
(see `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md` Phase 1c) - no new data
source, no change to how `bound_sd` itself is computed.

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
callers (and the Phase 1c kernel gating, which already treats missing
distinctness data as "no gating - fall back to the plain kernel") keep
working without OSD access.
