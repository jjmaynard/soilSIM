# Fit a distribution to a low/representative/high percentile triplet

The single dispatcher
[`extract_property_parameters()`](https://jjmaynard.github.io/soilSIM/reference/extract_property_parameters.md)
(`core-montecarlo.R`) calls to turn a SSURGO `_l/_r/_h` triplet into a
family-appropriate fit.

## Usage

``` r
fit_percentile_triplet(
  l,
  r,
  h,
  family,
  lh_probs = c(0.05, 0.95),
  bounds = NULL,
  boundedness = if (is.null(bounds)) "u" else "b"
)
```

## Arguments

- l, r, h:

  Low/representative/high values.

- family:

  One of `"triangular"`, `"uniform"`, `"normal"`, `"lognormal"`,
  `"beta"`, `"metalog"`, `"linear_cdf"`, `"auto"`.

- lh_probs:

  Length-2 vector giving the probabilities `l`/`h` are assumed to
  represent (SSURGO's low/high are professional-judgment bounds, not
  confirmed exact percentiles - this makes that assumption explicit and
  overridable rather than silently hardcoded).

- bounds:

  Optional `c(lower, upper)` physical bounds; required for `"beta"` and
  for `"metalog"` when `boundedness != "u"`.

- boundedness:

  Metalog boundedness; defaults to `"u"` when `bounds` is `NULL`, else
  `"b"`.

## Value

`list(family=, fit=, valid=)`. `valid=FALSE` when any of `l`/`r`/`h` is
non-finite.

## See also

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
[`compare_percentile_methods()`](https://jjmaynard.github.io/soilSIM/reference/compare_percentile_methods.md),
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md),
[`fit_beta_mle_newton()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton.md),
[`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md),
[`fit_normal_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_triplet.md),
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md),
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
[`quantile_beta()`](https://jjmaynard.github.io/soilSIM/reference/quantile_beta.md),
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
[`quantile_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf.md),
[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md),
[`quantile_triangular()`](https://jjmaynard.github.io/soilSIM/reference/quantile_triangular.md),
[`resolve_property_family()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_family.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
