# Recover SSURGO Components Missing All Horizon Data

[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)'s
`INNER JOIN chorizon` silently drops any real component that has zero
`chorizon` rows in SDA. This function finds such components (via
`all_components`, a component-only fetch with no chorizon join - see
[`fetch_ssurgo_all_components_working()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_all_components_working.md))
and, for each one, tries to synthesize a representative horizon profile
by averaging OTHER cokeys in the same AOI's `ssurgo_data` that share the
missing component's `compname` (see
[`synthesize_component_horizons_from_siblings()`](https://jjmaynard.github.io/soilSIM/reference/synthesize_component_horizons_from_siblings.md)).
Components with no such AOI sibling are left out of `ssurgo_data` (same
as current behavior) but are returned in `components_missing_horizons`
instead of vanishing untraced.

## Usage

``` r
recover_missing_horizon_components(
  ssurgo_data,
  all_components,
  verbose = FALSE
)
```

## Arguments

- ssurgo_data:

  The AOI's already fully processed horizon-grain data frame (i.e.
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
  Step 7 output - after RFV aggregation and restriction indicators,
  before
  [`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)).

- all_components:

  Data frame from
  [`fetch_ssurgo_all_components_working()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_all_components_working.md)
  (or shaped like it) - every real component for this AOI's
  `mukey_list`, independent of horizon data availability.

- verbose:

  Logical; provide progress messages.

## Value

`list(ssurgo_data = <ssurgo_data, with synthesized rows appended if any>, components_missing_horizons = <data frame: mukey/cokey/compname/comppct_l/r/h for components neither present in ssurgo_data nor resolvable via an AOI sibling>, components_recovered = <data frame: mukey/cokey/compname/n_sibling_cokeys_used for components that WERE synthesized>)`.
