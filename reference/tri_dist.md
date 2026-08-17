# Draw random samples from a triangular distribution

Ported verbatim from `code_ref/brdf/property_simulation.R`'s
`tri_dist()` (originally adapted from the `triangle` package), used by
`R/depth-simulation.R`'s profile-depth simulators. Kept as its own
random-draw implementation rather than layered on
[`quantile_triangular()`](https://jjmaynard.github.io/soilSIM/reference/quantile_triangular.md)
(a deterministic inverse-CDF evaluator with different degenerate-input
handling) so its exact edge-case behavior is preserved: it errors on
invalid `n`, but returns `NaN` (not an error) when the mode falls
outside `[a, b]` or any parameter is infinite/`NA`.

## Usage

``` r
tri_dist(n = 1, a = 0, b = 1, c = (a + b)/2)
```

## Arguments

- n:

  Number of samples to draw. If a vector is passed, its length is used
  (matching
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html)-style
  conventions).

- a, b, c:

  Minimum, maximum, and mode of the triangular distribution.

## Value

A numeric vector of length `n`; `NaN` for every element if `a`, `b`, `c`
describe a degenerate/invalid triangle.
