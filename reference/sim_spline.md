# Monotonic-spline inverse-CDF sampler (smoother than linear, still exact at knots)

Monotonic-spline inverse-CDF sampler (smoother than linear, still exact
at knots)

## Usage

``` r
sim_spline(probs, values, n, bounds = NULL)
```

## Arguments

- bounds:

  Optional length-2 vector giving the value at prob 0 and prob 1. If
  NULL, extends the observed value range by 5% on each side.
