# Method-of-moments Gamma fit

Used by the closed-form same-family fusion route's infeasible-cell
fallback (re-express as Normal moments, fuse as Normal, convert back).

## Usage

``` r
moments_to_gamma(mean, var)
```

## Arguments

- mean, var:

  Mean/variance to convert.

## Value

`list(shape=, rate=)`.
