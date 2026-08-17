# Rank fitted distributions by AIC (lower is better)

Rank fitted distributions by AIC (lower is better)

## Usage

``` r
rank_distribution_fits(fitted_distributions)
```

## Arguments

- fitted_distributions:

  Named list of fit results from fit_single_distribution() (NULL entries
  are dropped defensively, though fit_property_distributions() already
  filters them before this is called).

## Value

The same list, reordered best-fit-first.
