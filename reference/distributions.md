# Percentile-Triplet Distribution Fitting, ILR Composition Utilities, and Matrix Utilities

Shared foundation used by both `statistics.R` and `monte-carlo.R` (and,
for the ILR pieces, `bayesian-updating.R`).

The percentile-triplet fitting/quantile functions are closed-form base-R
vector arithmetic. The Beta fit is a vectorized Newton-Raphson MLE
(matches
[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
to ~1e-3 to 1e-5); the metalog fit is an exact linear solve when the
number of interior percentiles equals the number of terms (matches
`rmetalog::metalog()` to ~1e-12), so the `rmetalog` package is not
required.

The ILR (isometric log-ratio) transform matches `compositions::ilr()`
exactly, so the `compositions` package is not required either.
