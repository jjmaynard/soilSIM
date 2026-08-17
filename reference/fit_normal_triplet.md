# Fit a Normal distribution from three percentiles (exact closed form)

Fit a Normal distribution from three percentiles (exact closed form)

## Usage

``` r
fit_normal_triplet(p_lo_val, p50_val, p_hi_val, p_lo, p_hi)

quantile_normal(fit, q)
```

## Arguments

- p_lo_val, p50_val, p_hi_val:

  Values at probabilities `p_lo`, 0.5, `p_hi`.

- p_lo, p_hi:

  The low/high probabilities `p_lo_val`/`p_hi_val` represent.

- fit:

  A fit from `fit_normal_triplet()`.

- q:

  Vector of probabilities to evaluate the quantile function at.

## Value

`list(mean=, sd=)`.
