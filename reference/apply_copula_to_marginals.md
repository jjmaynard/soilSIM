# Map a Joint Copula Sample onto Existing Per-Depth Marginal Distributions

Converts the standard-normal joint sample from
[`sample_joint_depth_property_copula()`](https://jjmaynard.github.io/soilSIM/reference/sample_joint_depth_property_copula.md)
into actual property values, by probability-integral-transforming each
depth/property/realization cell
([`pnorm()`](https://rdrr.io/r/stats/Normal.html)) and then
inverse-transforming
([`quantile()`](https://rdrr.io/r/stats/quantile.html)) against that
depth/property's OWN already-simulated marginal distribution
(`property_matrices[[prop]][depth, ]` - the same per-horizon
triangular-fit values
[`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)
produced) - a Gaussian copula, in the standard statistical sense. This
is what actually delivers the achieved correlation structure from
[`sample_joint_depth_property_copula()`](https://jjmaynard.github.io/soilSIM/reference/sample_joint_depth_property_copula.md)
onto real property values while preserving each depth's own fitted
marginal shape exactly (for large `n_sims`).

## Usage

``` r
apply_copula_to_marginals(Z, property_matrices, gp_predictions = NULL)
```

## Arguments

- Z:

  A `c(n_depths, k, n_sims)` array as returned by
  [`sample_joint_depth_property_copula()`](https://jjmaynard.github.io/soilSIM/reference/sample_joint_depth_property_copula.md).

- property_matrices:

  Named list of `k` matrices (rows = depths, columns = simulations,
  `dim(Z)[1]`/`dim(Z)[3]` must match), in the same property order as
  `Z`'s second dimension - the same shape
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)
  already consumes.

- gp_predictions:

  Optional named list of per-property depth-trend mean vectors (length
  `n_depths`), e.g. from
  [`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md).
  When supplied for a property, that property's marginal at each depth
  is re-centered (a location shift only - shape/spread preserved) onto
  the GP-predicted mean before the quantile mapping, replacing the old
  sequential `gp_ratio` nudge with a single direct shift.

## Value

A named list of `n_depths x n_sims` matrices, same shape/names as
`property_matrices`.
