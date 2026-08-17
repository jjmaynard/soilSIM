# Estimate one horizon's ILR-space mean/covariance from marginal l/r/h triplets

Draws each fraction's own marginal Normal (clamped `> 0.01`), CLOSES
each draw (renormalizes the 3 draws to sum to 100 - the standard
compositional-data "closure" operation, which is what induces the
correct negative cross-fraction correlation, not an assumption),
ILR-transforms, and takes the empirical mean/covariance. A cheaper
closed-form delta-method shortcut was tried upstream and REJECTED: it
underestimates the true covariance by 5-10x once per-fraction CV exceeds
~30% (common near boundary fractions) - do not substitute it for this
Monte Carlo version.

## Usage

``` r
estimate_ilr_moments_mc(
  low_clay,
  rep_clay,
  high_clay,
  low_sand,
  rep_sand,
  high_sand,
  low_silt,
  rep_silt,
  high_silt,
  z,
  n_mc = 2000
)
```

## Arguments

- low_clay, rep_clay, high_clay, low_sand, rep_sand, high_sand,
  low_silt, rep_silt, high_silt:

  Low/representative/high values for each fraction.

- z:

  Standard-normal quantile matching the low/high interval (e.g.
  `qnorm(0.95)` for a 5th/95th-percentile low/high).

- n_mc:

  Monte Carlo sample size.

## Value

`list(mu = c(z1, z2), Sigma = 2x2 matrix)`.
