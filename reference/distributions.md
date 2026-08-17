# Percentile-Triplet Distribution Fitting, ILR Composition Utilities, and Matrix Utilities

Shared foundation used by both `statistics.R` and `monte-carlo.R` (and,
for the ILR pieces, `bayesian-updating.R`).

The percentile-triplet fitting/quantile functions below are adapted from
validated closed-form math originally written as
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
arithmetic
(`code_ref/reanalysis-platform/distribution_fitting_raster.R`): stripped
of all `terra::` calls, the underlying math is ordinary base-R vector
arithmetic, so it ports directly. The Beta fit (vectorized
Newton-Raphson MLE) was validated there against
[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
to ~1e-3 to 1e-5; the metalog fit (exact linear solve when the number of
interior percentiles equals the number of terms) was validated to ~1e-12
against `rmetalog::metalog()` - this project therefore never needs the
`rmetalog` package (which `code_ref/brdf/distribution_fitting.R`
documents as prone to hangs and a reproduced segfault).

The ILR (isometric log-ratio) transform is ported verbatim from
`code_ref/reanalysis-platform/texture_ilr_fusion.R`, itself validated to
match `compositions::ilr()` exactly - so this package needs no
`compositions` dependency either.
