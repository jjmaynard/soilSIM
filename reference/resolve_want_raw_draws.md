# Resolve whether a `run_stage1_fusion()` call should fuse against real per-mukey Monte Carlo draws (`prior_fusion_method = "raw_draws"`) or the percentile-reconstructed approximation.

`MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` task P2.10: scope-aware default
flip. An explicit `prior_fusion_method` value ("raw_draws" or
"percentile") always wins and is honored exactly as before - this
function only changes what happens when it is unset (`NULL`, the common
case).

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

Unset now defaults to raw_draws for any `dist` whose fusion route
actually has a raw-draws fit implemented: `"normal"`/`"beta"`/`"gamma"`
(both AOI-size routes, as of P2.6) and `"lognormal"` (both its small-AOI
general branch and large-AOI closed-form branch, as of P2.6/P2.12).
`"auto"` is included too:
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md)'s
auto-selection only ever resolves to `"beta"`/`"normal"`/`"lognormal"`
(see its own code) - **never** `"metalog"` - so it's safe to default in
as well. Only an explicit `dist = "metalog"` keeps the original
percentile- reconstruction default - **not** because no raw-draws fit
exists for it (P2.12 added one, an empirical-percentile variant of
[`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)'s
own interpolation), but because its cost/benefit hasn't been
benchmarked, exactly mirroring the texture-group route's P2.11
treatment: that default-flip decision is deliberately kept separate from
"does a fit exist."

This default was set only after confirming the added cost is modest at
any AOI size: raw_draws costs 16-25% more than the default on the
general-KDE route (P2.9) and ~0.18s of overhead on top of the
closed-form route's own ~19.73s fit cost at 100,172 cells/150 mukeys
(P2.6/P2.10 investigation) - i.e. this flip trades a small, bounded cost
increase for fusing against real simulated data instead of a 5-9 point
percentile reconstruction, for every property that can use it, without
the caller having to know to ask for it.
