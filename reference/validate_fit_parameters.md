# Validate a fitted distribution's parameters

Real replacement for the always-`valid=TRUE`
[`validate_distribution_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_distribution_parameters.md)
stub previously in `mod05_monte_carlo.R`.

## Usage

``` r
validate_fit_parameters(family, fit)
```

## Arguments

- family:

  One of
  [`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)'s
  resolved families.

- fit:

  The `fit` element of
  [`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)'s
  return value.

## Value

`list(valid=, message=)`.
