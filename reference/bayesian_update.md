# Combine a prior and likelihood distribution into a posterior via grid-based KDE

Estimates prior and likelihood densities over a shared value grid via
kernel density estimation, multiplies them per Bayes' rule, and samples
from the resulting (normalized) posterior. Fully general - no
distributional family assumption, no requirement that both sides match -
at the cost of needing raw samples (not distributional parameters) and
being considerably more expensive than the closed-form tiers above.

## Usage

``` r
bayesian_update(
  prior_distribution,
  likelihood_distribution,
  grid_range = NULL,
  grid_resolution = 0.01,
  n = 1000
)
```

## Arguments

- prior_distribution, likelihood_distribution:

  Numeric vectors of samples.

- grid_range:

  Optional length-2 vector giving the grid's min/max. If `NULL`, derived
  from the combined range of both inputs, padded by 1.

- grid_resolution:

  Step size of the evaluation grid.

- n:

  Number of posterior samples to draw.

## Value

A numeric vector of `n` samples drawn from the posterior distribution.
