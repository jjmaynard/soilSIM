# Closed-form same-family fusion path for `fuse_adaptive()`, with a per-cell fallback to Normal-moment fusion (re-expressed back into the requested family) where the naive same-family fusion is infeasible.

Closed-form same-family fusion path for
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
with a per-cell fallback to Normal-moment fusion (re-expressed back into
the requested family) where the naive same-family fusion is infeasible.

## Usage

``` r
fuse_closed_form(
  prior_value_rasters,
  lik_value_rasters,
  percentile_probs,
  family,
  bounds
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Lists of percentile-value SpatRasters.

- percentile_probs:

  Numeric probabilities matching the value-raster list order.

- family:

  One of "normal", "beta", "gamma".

- bounds:

  Required for `family = "beta"`.

## Value

`list(posterior=, route_detail=, n_fallback_cells=)`.
