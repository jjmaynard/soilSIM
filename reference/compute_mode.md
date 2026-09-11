# Compute the Mode of a Numeric Vector

Compute the Mode of a Numeric Vector

## Usage

``` r
compute_mode(x)
```

## Arguments

- x:

  A vector of values.

## Value

The most frequently-occurring value in `x`.

## Behavior

`tabulate(match(x, ux))` on pre-sorted unique values counts occurrences
without the factor-coercion overhead of `table(x)`/`factor(x)`. On ties
it returns the smallest value among them (unique values are sorted
ascending and [`which.max()`](https://rdrr.io/r/base/which.min.html)
returns the first maximum).
