# Evaluate a metalog quantile function

Unlike the raster-native original (which evaluates one shared `q` across
many cells' coefficient rasters), here `fit` is a single set of
coefficients and `q` may be a vector - so this multiplies the basis
matrix evaluated at `q` by the coefficient vector directly, rather than
the raster version's "select row 1" pattern.

## Usage

``` r
quantile_metalog_linear(fit, q)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md).

- q:

  Vector of probabilities in (0,1).
