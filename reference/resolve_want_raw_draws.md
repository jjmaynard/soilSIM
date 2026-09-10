# Resolve whether a `run_stage1_fusion()` call should fuse against real per-mukey Monte Carlo draws (`prior_fusion_method = "raw_draws"`) or the percentile-reconstructed approximation.

An explicit `prior_fusion_method` value ("raw_draws" or "percentile")
always wins; this function only decides what happens when it is unset
(`NULL`, the common case).

## Usage

``` r
resolve_want_raw_draws(prior_fusion_method, dist)
```

## Arguments

- prior_fusion_method:

  `property_config$prior_fusion_method` - `NULL`/unset resolves per
  `dist` (see above); `"raw_draws"`/`"percentile"` are explicit
  overrides, always honored.

- dist:

  `property_config$dist` - the value as configured, BEFORE
  [`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md)
  resolves `"auto"` to a concrete family (this function runs earlier,
  before the mukey draws simulation that resolution would otherwise
  duplicate the cost of triggering).

## Value

`TRUE`/`FALSE`.

## Details

Unset defaults to raw_draws for any `dist` whose fusion route has a
raw-draws fit: `"normal"`/`"beta"`/`"gamma"` (both AOI-size routes) and
`"lognormal"` (both its small-AOI general branch and large-AOI
closed-form branch). `"auto"` is included, since
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md)'s
auto-selection only ever resolves to `"beta"`/`"normal"`/ `"lognormal"`,
never `"metalog"`. An explicit `dist = "metalog"` uses the percentile
reconstruction:
[`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)
has an empirical-percentile raw-draws variant, but its cost/benefit is
not benchmarked, so it is kept off by default.

On the general-KDE route raw_draws costs about 16-25% more than the
percentile path; on the closed-form route it adds about 0.18s on top of
the fit cost at 100,172 cells / 150 mukeys.
