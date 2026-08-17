# Vectorized Newton-Raphson Beta MLE from raw-scale percentile values

Fits a single `(shape1, shape2)` pair from `k` percentile values of one
property at one horizon. The Beta log-likelihood's score/Hessian have a
closed form in
[`digamma()`](https://rdrr.io/r/base/Special.html)/[`trigamma()`](https://rdrr.io/r/base/Special.html)
(standard Beta MLE theory), iterated for a fixed number of steps rather
than via a general-purpose optimizer - validated (upstream) against
[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
to ~1e-3-1e-5 after 15 iterations.

`fit_beta_mle_newton_vec()` runs the identical update across every row
of `value_matrix` simultaneously
([`digamma()`](https://rdrr.io/r/base/Special.html)/[`trigamma()`](https://rdrr.io/r/base/Special.html)
are already vectorized base-R functions), so it is the version
[`prepare_simulation_parameters()`](https://jjmaynard.github.io/soilSIM/reference/prepare_simulation_parameters.md)
should call once per property across all horizons, rather than looping
`fit_beta_mle_newton()` per row.

## Usage

``` r
fit_beta_mle_newton(values, bounds, n_iter = 15, eps = 1e-06)

fit_beta_mle_newton_vec(value_matrix, bounds, n_iter = 15, eps = 1e-06)
```

## Arguments

- values:

  Numeric vector of `k` percentile values (raw units).

- bounds:

  `c(lower, upper)` physical bounds.

- n_iter:

  Fixed Newton-Raphson iteration count.

- eps:

  Clamp distance from the exact 0/1 rescaled boundary (avoids `-Inf` in
  the score equations when a percentile lands exactly on a bound).

- value_matrix:

  Numeric matrix, one row per horizon, one column per percentile value
  (e.g. 3 columns for an l/r/h triplet).

## Value

`list(shape1=, shape2=)`.

A data frame with columns `shape1`, `shape2` (one row per horizon).
