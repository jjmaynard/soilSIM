# Resolve a property's distribution family from its own percentile skew

Computed directly from this row's own l/r/h (no spatial/AOI aggregation
in a tabular context). Deliberately narrow - never resolves to
`"metalog"` (that stays an explicit config choice), matching the
reference's own documented restriction.

## Usage

``` r
resolve_property_family(l, r, h, bounds = NULL, skew_threshold = 0.15)
```

## Arguments

- l, r, h:

  Low/representative/high values.

- bounds:

  Optional `c(lower, upper)`; when supplied, always resolves to
  `"beta"`.

- skew_threshold:

  Absolute skew-proxy threshold below which `"normal"` is chosen over
  `"lognormal"`.

## Value

`list(family=, skew_proxy=)`.

## See also

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
[`compare_percentile_methods()`](https://jjmaynard.github.io/soilSIM/reference/compare_percentile_methods.md),
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md),
[`fit_beta_mle_newton()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton.md),
[`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md),
[`fit_normal_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_triplet.md),
[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md),
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md),
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
[`quantile_beta()`](https://jjmaynard.github.io/soilSIM/reference/quantile_beta.md),
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
[`quantile_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf.md),
[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md),
[`quantile_triangular()`](https://jjmaynard.github.io/soilSIM/reference/quantile_triangular.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
