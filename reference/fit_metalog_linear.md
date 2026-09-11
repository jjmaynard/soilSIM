# Fit a metalog distribution via exact linear solve

When the number of interior percentiles (`p` strictly in `(0,1)`) equals
the number of metalog terms, the fit is an EXACTLY-DETERMINED linear
system (`solve(Y, z)`), not an optimization - validated (upstream) to
~1e-12 against `rmetalog::metalog()`. This is the fast path `rmetalog`
itself takes when its solution is already feasible; see
[`check_metalog_feasible()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasible.md)/[`quantile_metalog_linear_fallback()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_fallback.md)
for the guard `rmetalog`'s own LP feasibility-correction would otherwise
provide.

## Usage

``` r
fit_metalog_linear(
  interior_values,
  interior_probs,
  bounds = NULL,
  boundedness = "u"
)
```

## Arguments

- interior_values, interior_probs:

  Percentile values/probabilities, both strictly interior (excluding
  p=0/p=1 if present).

- bounds:

  `c(lower, upper)`; required unless `boundedness="u"`.

- boundedness:

  One of `"u"` (unbounded), `"sl"` (semi-bounded below), `"su"`
  (semi-bounded above), `"b"` (bounded).

## Value

`list(a=, term=, bounds=, boundedness=)`.

## See also

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
[`compare_percentile_methods()`](https://jjmaynard.github.io/soilSIM/reference/compare_percentile_methods.md),
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md),
[`fit_beta_mle_newton()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton.md),
[`fit_normal_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_triplet.md),
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
