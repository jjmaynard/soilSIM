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

## See also

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
[`compare_percentile_methods()`](https://jjmaynard.github.io/soilSIM/reference/compare_percentile_methods.md),
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md),
[`fit_beta_mle_newton()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton.md),
[`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md),
[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md),
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md),
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
[`quantile_beta()`](https://jjmaynard.github.io/soilSIM/reference/quantile_beta.md),
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
[`quantile_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf.md),
[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md),
[`quantile_triangular()`](https://jjmaynard.github.io/soilSIM/reference/quantile_triangular.md),
[`resolve_property_family()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_family.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
[`compare_percentile_methods()`](https://jjmaynard.github.io/soilSIM/reference/compare_percentile_methods.md),
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md),
[`fit_beta_mle_newton()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton.md),
[`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md),
[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md),
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md),
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
[`quantile_beta()`](https://jjmaynard.github.io/soilSIM/reference/quantile_beta.md),
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
[`quantile_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf.md),
[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md),
[`quantile_triangular()`](https://jjmaynard.github.io/soilSIM/reference/quantile_triangular.md),
[`resolve_property_family()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_family.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
