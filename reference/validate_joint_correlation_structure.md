# Validate Joint Depth x Property Correlation Structure

Diagnostic helper for the joint depth x property correlation structure.
Unlike
[`validate_correlation_preservation()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_preservation.md),
which only checks cross-property correlation and (via its caller
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md))
only at the first 5 depths, this reports BOTH halves of the joint
structure - cross-property correlation at every depth, and depth-lag
(across-depth, within-property) correlation for every property - against
explicit target matrices when supplied. This is the acceptance test any
vertical-correlation method (the existing GP-quantile-retrofit path or a
future joint-copula replacement) should be judged against: does the SAME
simulated output actually achieve both the intended property correlation
matrix and the intended depth correlation structure, simultaneously.

## Usage

``` r
validate_joint_correlation_structure(
  simulated_list,
  target_property_corr = NULL,
  target_depth_corr = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulated_list:

  A named list of matrices (rows = depths, columns = simulations), the
  same shape consumed/produced by
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)/
  [`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md).

- target_property_corr:

  Optional k x k target property-correlation matrix (row/column names
  matching `names(simulated_list)`, or unnamed matching order). If
  supplied, per-depth and overall max absolute deviation from this
  target is reported.

- target_depth_corr:

  Optional n_depths x n_depths target depth-correlation matrix. If
  supplied, per-property and overall max absolute deviation from this
  target is reported.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE`

  - quiet). See
    [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

A list with elements:

- achieved_property_correlation:

  List (one per depth) of achieved k x k property correlation matrices.

- property_correlation_max_diff:

  Numeric vector (one per depth) of max absolute deviation from
  `target_property_corr`, or all `NA` if no target supplied.

- achieved_depth_correlation:

  Named list (one per property) of achieved n_depths x n_depths depth
  correlation matrices.

- depth_correlation_max_diff:

  Named numeric vector (one per property) of max absolute deviation from
  `target_depth_corr`, or all `NA` if no target supplied.

- overall_property_max_diff, overall_depth_max_diff:

  Single worst-case summary values across all depths/properties
  respectively.
