# Raster-Native Percentile-Based Distribution Fitting

Fits a distribution to a property's percentile-value rasters (e.g.
SSURGO low/rep/high, SOLUS low/pred/high) as pure
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
arithmetic and/or
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html)/[`terra::lapp()`](https://rspatial.github.io/terra/reference/lapp.html)
calls to *vectorized* base-R math functions (`qnorm`, `qbeta`,
`digamma`/`trigamma`) - never a per-cell optimizer, which is what makes
these fast across whole rasters. Ported from
`code_ref/reanalysis-platform/distribution_fitting_raster.R`.

Reuses `R/distributions.R`'s existing
[`metalog_basis_matrix()`](https://jjmaynard.github.io/soilSIM/reference/metalog_basis_matrix.md)/
`metalog_to_z()`/`metalog_from_z()` directly rather than duplicating
them: those functions are pure elementwise arithmetic (`log`/`exp`/`/`
plus building a small, non-spatial coefficient matrix from the fixed
probability grid), so they already work unchanged whether their `x`/`z`
arguments are plain numerics or
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)s -
confirmed by a dedicated smoke test before this file was written, not
assumed.

Status, from the original source:

- [`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)
  /
  [`quantile_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md) -
  EXACT.

- [`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md) -
  EXACT (nonparametric fallback target).

- [`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md) -
  method-of-moments Beta fit, APPROXIMATE relative to MLE, used only as
  the Newton-Raphson seed below.

- [`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md)
  /
  [`quantile_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md) -
  vectorized Newton-Raphson Beta MLE, validated (upstream) directly
  against
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
  (~1e-3 to 1e-5).

- [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
  /
  [`quantile_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
  /
  [`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)
  /
  [`quantile_metalog_linear_with_fallback()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback.md) -
  exact linear-solve metalog reformulation (matches
  `rmetalog::metalog()` to ~1e-12 when its own fit is feasible, per
  upstream validation), with a validated fallback to `linear_cdf` for
  the cells where it isn't.

`spline` is intentionally excluded (per the original project's own
decision - its accuracy relative to R's `splinefun(method="monoH.FC")`
was never confirmed bit-exact).
