# Classify a horizon name into a generalized master-horizon designation

Strips a leading lithologic-discontinuity digit prefix (e.g. `"2Bt2"`
-\> `"Bt2"`), then matches against the master-horizon letters the KSSL
reference correlation matrices are keyed by: `O`, `A`, `E`, `B`, `C`,
`Cr`, `R`. `Cr` is checked before `C` so `"Cr"`/`"Crt"` aren't
misclassified as `"C"`.

## Usage

``` r
classify_genhz(hzname)
```

## Arguments

- hzname:

  Character vector of raw SSURGO horizon-designation text (e.g. `"Bt2"`,
  `"2Bw"`, `"Cr"`, `"R"`). `NA` is preserved as `NA`.

## Value

Character vector, same length as `hzname`, of
`O`/`A`/`E`/`B`/`C`/`Cr`/`R` or `NA_character_` for unrecognized text.

## See also

Other distributions:
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
[`resolve_property_family()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_family.md),
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
