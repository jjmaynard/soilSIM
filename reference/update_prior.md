# Combine a prior and likelihood distribution into a posterior via grid-based KDE

Estimates prior and likelihood densities over a shared value grid via
kernel density estimation, multiplies them per Bayes' rule, and samples
from the resulting (normalized) posterior. Fully general - no
distributional family assumption, no requirement that both sides match -
at the cost of needing raw samples (not distributional parameters) and
being considerably more expensive than the closed-form tiers above.

## Usage

``` r
update_prior(
  prior_distribution,
  likelihood_distribution,
  grid_range = NULL,
  grid_resolution = 0.01,
  n = 1000,
  winsorize_probs = NULL,
  posterior_probs = NULL
)
```

## Arguments

- prior_distribution, likelihood_distribution:

  Numeric vectors of samples.

- grid_range:

  Optional length-2 vector giving the grid's min/max. If `NULL`, derived
  from the combined range of both inputs, padded by 1.

- grid_resolution:

  Step size of the evaluation grid.

- n:

  Number of posterior samples to draw.

- winsorize_probs:

  Optional length-2 vector `c(lower, upper)` (e.g. `c(0.01, 0.99)`).
  When supplied, each side's sample is independently clipped to its OWN
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) range at these
  probabilities before the `stats::density(bw = "nrd0")` call. `NULL`
  (default) preserves this function's original behavior exactly.

- posterior_probs:

  Optional numeric vector of probabilities (0-1) at which to report
  posterior percentiles (e.g. `c(0.05, 0.5, 0.95)`). `NULL` (default)
  preserves this function's original return shape EXACTLY - a plain
  numeric vector of `n` posterior samples - so every existing caller
  ([`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
  `R/core-montecarlo.R`'s
  [`fuse_observed_data_into_priors()`](https://jjmaynard.github.io/soilSIM/reference/fuse_observed_data_into_priors.md),
  this file's own tests) is completely unaffected. When supplied, the
  return shape changes to the richer list described below (see
  `@return`); see `@section Grid-based percentiles` for why this is the
  primary use case for `posterior_probs`, not
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) on the resampled
  vector.

## Value

When `posterior_probs` is `NULL` (default): a numeric vector of `n`
samples drawn from the posterior distribution (unchanged from this
function's original contract). When `posterior_probs` is supplied: a
list with -

- `samples`: the same `n`-length resampled vector as the `NULL` case
  (still available for callers that want actual posterior draws, e.g.
  downstream Monte Carlo propagation).

- `mean`, `var`: the posterior mean/variance computed EXACTLY from the
  discretized `posterior_prob`/`value_grid` distribution
  (`sum(value_grid * posterior_prob)` and its variance analogue) - zero
  resampling noise, strictly at least as accurate as `mean(samples)`/
  `var(samples)` for any `n`.

- `percentiles`: a named numeric vector (names `"P5"`, `"P50"`, ...
  matching `posterior_probs`), read directly off
  `cumsum(posterior_prob)` - also zero resampling noise, accurate to
  within `grid_resolution`.

- `value_grid`, `posterior_prob`: the underlying discretized
  distribution itself, for callers that want to compute something these
  three summaries don't cover.

## Grid-based percentiles have zero resampling noise

This function already computes the FULL discretized posterior
distribution (`posterior_prob` over `value_grid`) before its final step,
which historically only existed to hand back a plain vector
(`sample(value_grid, size = n, prob = posterior_prob, replace = TRUE)`).
Reading percentiles/moments directly off `posterior_prob`'s cumulative
distribution - finding where `cumsum(posterior_prob)` crosses each
target probability - is exact relative to `grid_resolution` and has no
resampling noise at all, unlike `quantile(samples, probs = ...)` would
have. No amount of increasing `n` improves on reading percentiles
directly off an already-fully-known discrete distribution, so this is
`posterior_probs`' primary computation path, not a convenience wrapper
around the resampled `samples` vector.

## Known limitation (raw-draws fusion vs. percentile-reconstruction fusion)

A dedicated adversarial test found that for a skewed prior fused against
a weak/wide likelihood, `soilSIM`'s `prior_fusion_method = "raw_draws"`
route (a genuine empirical resample as the prior side) can produce a
fused posterior mean measurably *less* close to the prior population's
true mean than the default percentile-reconstruction route, in a way
that first looked like a bandwidth-selection bug. Task P3.1's follow-up
investigation ruled that out directly (forcing both routes onto the
exact same [`stats::density()`](https://rdrr.io/r/stats/density.html)
bandwidth left the gap essentially unchanged) and traced the real
mechanism to `R/core-distributions.R`'s
[`sim_linear_cdf_batch()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf_batch.md):
reconstructing a distribution from only 5 percentile knots is an
inherent information-loss approximation for skewed shapes (confirmed
method-agnostic - `linear_cdf`/`spline`/`kde` reconstruction all show a
comparable P75-P95 overshoot on the same synthetic scenario), not a
defect in this function or in raw_draws. See that task's write-up for
the full investigation; `winsorize_probs` remains available below as a
harmless, independently-useful safety net for genuine extreme-outlier
inputs, not as a fix for this finding (there was nothing here to fix).

## See also

Other bayesian-fusion:
[`convert_beta_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_beta_to_moments.md),
[`convert_gamma_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_gamma_to_moments.md),
[`convert_lognormal_to_normal()`](https://jjmaynard.github.io/soilSIM/reference/convert_lognormal_to_normal.md),
[`convert_moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_beta.md),
[`convert_moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_gamma.md),
[`convert_normal_to_lognormal()`](https://jjmaynard.github.io/soilSIM/reference/convert_normal_to_lognormal.md),
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md),
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md),
[`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md),
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md)
