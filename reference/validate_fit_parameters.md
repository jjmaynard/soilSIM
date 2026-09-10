# Validate a fitted distribution's parameters

Checks that a fitted distribution's parameters are structurally complete
and in range for its family.

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
