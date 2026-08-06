# Fit a Single Local GP Model

Fit a Single Local GP Model

## Usage

``` r
fit_local_gp_model_single(agg_data, prop, gp_control = c(20, 10, 2))
```

## Arguments

- agg_data:

  One row per depth (`hzdept_r`, `mean_val`), as returned by
  `aggregate_property_by_depth()` - always a handful of points
  (typically \<10, one per unique depth in a single cokey), never one
  row per Monte Carlo replicate.

- prop:

  Property name, stored on the returned model for reference.

- gp_control:

  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s
  `control` argument (population size / iteration counts for its
  internal hyperparameter search). Defaults to `c(20, 10, 2)`, far below
  `GP_fit()`'s own default `c(200*d, 80*d, 2*d)` (`d` = input
  dimensionality, always 1 here - depth is the only predictor).
  `GP_fit()`'s default search effort scales with `d`, not with the
  number of training points, and a 1-D, single-hyperparameter (`beta`)
  fit on this few points has a simple enough likelihood surface that the
  smaller search converges to the identical optimum every time -
  verified empirically across several representative depth/value series
  (identical fitted `beta` and identical predictions vs. the default,
  ~5x faster per call). Pass `c(200, 80, 2)` (or larger) to restore
  `GP_fit()`'s own default search effort if a future property/dataset
  needs a more thorough search.

## Value

A list with the fitted GP model, depth scaling info, training data, and
`prop`.
