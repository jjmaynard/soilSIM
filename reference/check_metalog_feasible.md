# Check whether a metalog fit's quantile function is monotonic (feasible)

Probes
[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md)
at a grid of probabilities and flags non-monotonicity (equivalent to the
implied density going negative somewhere) - a probe, not a proof, but
cheap and effective.

## Usage

``` r
check_metalog_feasible(fit, y_grid = seq(0.02, 0.98, by = 0.02))
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md).

- y_grid:

  Probability grid to probe.

## Value

`TRUE` if infeasible (non-monotonic), else `FALSE`.
