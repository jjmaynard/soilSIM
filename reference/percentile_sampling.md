# Flexible Percentile-Based Distribution Sampling

Consolidates several "simulate a distribution from summary percentiles"
approaches (a piecewise-linear inverse-CDF, a monotonic-spline
inverse-CDF, a KDE/truncated-normal approach, and beta/normal parametric
fits) behind a single `method` argument, for arbitrary percentile counts
(not just SSURGO's low/rep/high triplet). Ported from
`code_ref/brdf/distribution_fitting.R`. Complements, rather than
replaces, `R/distributions.R`'s
[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)/
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
which is the production Monte Carlo engine's fixed 3-point fitter.

`method = "metalog"` is deliberately **not** ported: it would need the
`rmetalog` package, which `distributions.R`'s own header documents a
project decision to avoid entirely (validated hang/segfault history) -
soilSIM already has a dependency-free closed-form metalog fitter
([`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md)/[`quantile_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear.md))
that doesn't need it.
