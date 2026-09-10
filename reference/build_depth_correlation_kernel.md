# Build a Depth Correlation Kernel Matrix

Constructs an `n_depths x n_depths` correlation matrix from real depth
*distances* (not just adjacent-lag steps), via an exponential or Matern
kernel with a given `length_scale`. This is the depth half of a
Kronecker-separable `R_depth (x) R_property` joint covariance (see
[`sample_joint_depth_property_copula()`](https://jjmaynard.github.io/soilSIM/reference/sample_joint_depth_property_copula.md)).

## Usage

``` r
build_depth_correlation_kernel(
  depths,
  length_scale,
  kernel = c("exponential", "matern"),
  nu = 1.5,
  boundary_distinctness = NULL,
  distinctness_range = c(min = 1, max = 10),
  min_gate_weight = 0.05
)
```

## Arguments

- depths:

  Numeric vector of depth values (real units, e.g. cm; need not be
  evenly spaced or sorted - internally sorted for the
  `boundary_distinctness` gating pass, then mapped back to `depths`'
  original order for the returned matrix).

- length_scale:

  A single positive numeric length-scale (real units matching `depths`),
  typically from
  [`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md).
  A non-finite or non-positive value degrades gracefully to the identity
  matrix (no depth correlation) rather than erroring, matching this
  package's established fallback pattern for missing/invalid correlation
  inputs (see
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
  pooled-matrix fallback).

- kernel:

  `"exponential"` (default) or `"matern"`.

- nu:

  Matern smoothness parameter (only used when `kernel = "matern"`);
  default `1.5`.

- boundary_distinctness:

  Optional numeric vector, same length as `depths`, of OSD-derived
  `bound_sd` values (see
  [`attach_osd_boundary_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/attach_osd_boundary_distinctness.md)) -
  one value per depth, interpreted as the distinctness of the boundary
  immediately ABOVE that depth (in sorted-ascending order; the value for
  the shallowest depth is unused, since there is no boundary above the
  top of the profile). Smaller `bound_sd` (e.g.
  [`aqp::hzDistinctnessCodeToOffset()`](https://ncss-tech.github.io/aqp/reference/hzDistinctnessCodeToOffset.html)'s
  `abrupt` ~= 1, `clear` ~= 2.5) means a sharper real discontinuity and
  more strongly suppresses correlation across that boundary; larger
  `bound_sd` (`gradual` ~= 7.5, `diffuse` ~= 10) leaves the plain
  distance-decay kernel almost untouched. `NULL` (default) or all-`NA`
  skips gating entirely, reproducing the plain distance-decay kernel
  exactly.

- distinctness_range:

  Named `c(min=, max=)` giving the `bound_sd` values that map to the
  strongest (`min_gate_weight`) and weakest (`1`, no suppression) gating
  respectively. Default `c(min=1, max=10)` spans
  [`aqp::hzDistinctnessCodeToOffset()`](https://ncss-tech.github.io/aqp/reference/hzDistinctnessCodeToOffset.html)'s
  actual `abrupt`-to-`diffuse` range for the codes
  [`infill_missing_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_distinctness.md)
  ever assigns.

- min_gate_weight:

  The gating floor (default `0.05`) - even the sharpest boundary leaves
  a small residual correlation rather than exactly zero, both physically
  (a truly zero-correlation entry is rarely justified) and numerically
  (helps
  [`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md)'s
  eigenvalue repair stay well-conditioned).

## Value

A symmetric, unit-diagonal, positive-definite
`length(depths) x length(depths)` correlation matrix (repaired via
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md)
if numerical floating-point asymmetry from the distance-matrix/gating
construction pushes it slightly off).
