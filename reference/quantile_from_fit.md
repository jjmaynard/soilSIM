# Evaluate the quantile function of a fitted percentile-triplet distribution

Replaces the old `transform_to_distribution()` design: a thin dispatcher
from `family` to the matching `quantile_*()` function.

## Usage

``` r
quantile_from_fit(u, family, fit)
```

## Arguments

- u:

  Vector of probabilities (e.g. correlated uniform draws).

- family:

  One of
  [`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)'s
  resolved families.

- fit:

  The `fit` element of
  [`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)'s
  return value.

## Value

Numeric vector, same length as `u`. `NA_real_` (not a silent wrong value
or error) when `fit` is `NULL`/invalid.
