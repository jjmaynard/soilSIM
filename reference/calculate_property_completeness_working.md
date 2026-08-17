# Calculate Property Completeness

Fraction of non-`NA` `_r` values per known soil property in `df`.

## Usage

``` r
calculate_property_completeness_working(df)
```

## Arguments

- df:

  Horizon data frame.

## Value

Named list of `property -> completeness fraction` (`[0, 1]`).
