# Default posterior percentile set computed by `fuse_general_kde()` (and, as later routes adopt it, every other fusion route).

Default posterior percentile set computed by
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
(and, as later routes adopt it, every other fusion route).

## Usage

``` r
FUSE_POSTERIOR_DEFAULT_PROBS
```

## Format

An object of class `numeric` of length 9.

## Why 9 points, not the 5-point SSURGO/SOLUS input convention

The 5-point set (`c(0.05, 0.25, 0.5, 0.75, 0.95)`) is what SSURGO/SOLUS
happen to provide as raw *input* percentile data - a real external
constraint on the prior/likelihood side. The *posterior* has no such
constraint, since this package computes it, and 5 points under-resolves
exactly the tail behavior risk/cost-assessment use cases need
(exceedance probabilities, asymmetric credible intervals). This set
covers a 98% interval (P01/P99, standard for tail-risk work), a 90%
interval (P05/P95), and an 80% interval (P10/P90), alongside the
existing median/IQR.

## Reliability note

For the closed-form analytic routes (`qnorm`/`qbeta`/`qgamma`), extra
percentiles are free and exact regardless of how extreme they are. For
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
sample/grid-based route, percentiles are read directly off
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
discretized `posterior_prob` (see that function's own "Grid-based
percentiles" section) - exact relative to `grid_resolution` and NOT
resampling-noise-limited, but still subject to how well the underlying
KDE itself approximates the tails from a finite `n_samples`-sized input
resample; the extreme percentiles here (P01/P99) are inherently noisier
than the inner ones (P10-P90) for that reason, not because of
resampling.

A caller needing even more extreme tail resolution (e.g. P0.1/P99.9) can
pass its own `posterior_probs` instead of relying on this default.
