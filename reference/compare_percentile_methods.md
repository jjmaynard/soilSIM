# Run several percentile-reconstruction methods on the same data and compare them

Draws `n` samples from each requested method and returns both the raw
samples and a summary-statistics table (via
[`compute_summary_statistics()`](https://jjmaynard.github.io/soilSIM/reference/compute_summary_statistics.md))
side by side, so the shape/spread of each reconstruction can be compared
by eye or by downstream tests (e.g.
[`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html) between pairs
of methods).

## Usage

``` r
compare_percentile_methods(
  quantile_df,
  methods = c("linear_cdf", "spline", "kde"),
  percentile_cols = c("P0", "P5", "P50", "P95", "P100"),
  n = 1000,
  bounds = NULL,
  ...
)
```

## Arguments

- quantile_df:

  A data frame with at least one row and percentile-named columns (e.g.
  "P0", "P5", "P50", "P95", "P100").

- methods:

  Character vector of methods to run (see `simulate_from_percentiles`).

- percentile_cols:

  Character vector of candidate percentile column names.

- n:

  Number of samples to draw.

- bounds:

  Optional length-2 vector of (lower, upper) bounds. Required for method
  = "beta"; used as padding bounds for "spline" if given.

- ...:

  Additional method-specific arguments (`sample_size` for "kde").

## Value

A list with `samples` (named list of numeric vectors, one per method)
and `summary` (data frame of summary statistics, one row per method).

## See also

Other distributions:
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md),
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
[`resolve_property_family()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_family.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
