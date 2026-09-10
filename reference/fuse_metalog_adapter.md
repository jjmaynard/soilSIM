# Metalog adapter for `fuse_property_adaptive()`: closed-form route via metalog-fitted Normal-equivalent moments, used regardless of AOI size.

Metalog adapter for
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md):
closed-form route via metalog-fitted Normal-equivalent moments, used
regardless of AOI size.

## Usage

``` r
fuse_metalog_adapter(
  prior_value_rasters,
  prior_probs,
  lik_value_rasters,
  lik_probs,
  property_config,
  verbose = TRUE,
  posterior_probs = NULL,
  mukey_raster = NULL,
  mukey_draws = NULL,
  ...
)
```

## Arguments

- posterior_probs:

  Numeric probabilities (0-1) at which to report posterior percentiles.
  `NULL` (default) resolves to
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md).
  Computed via
  [`qnorm_percentiles_raster()`](https://jjmaynard.github.io/soilSIM/reference/qnorm_percentiles_raster.md)
  on the final fused `(mu, sigma)` - this route already reduces both
  sides to Normal-equivalent moments
  ([`metalog_moments_raster()`](https://jjmaynard.github.io/soilSIM/reference/metalog_moments_raster.md)'s
  quadrature) before the final
  [`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
  step, so there is no richer distributional shape to draw percentiles
  from beyond that Normal approximation; carries forward this route's
  existing "less-validated glue code" caveat (see
  [`metalog_moments_raster()`](https://jjmaynard.github.io/soilSIM/reference/metalog_moments_raster.md)'s
  own docs) to its percentile output too.

- mukey_raster, mukey_draws:

  Optional raw-draws inputs. Metalog has no moment/density fit analogous
  to
  [`fuse_closed_form()`](https://jjmaynard.github.io/soilSIM/reference/fuse_closed_form.md)'s
  normal/beta/gamma/lognormal routes (its fit is an exact linear
  interpolation through fixed percentile knots) - the real extension is
  interpolating through each mukey's REAL empirical percentiles
  ([`mukey_draws_percentiles_raster()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_percentiles_raster.md))
  instead of the percentile-reconstruction values, wherever a mukey has
  draws coverage. `NULL` (default, or absorbed via `...` if passed
  positionally-incompatible) preserves the original
  percentile-reconstruction behavior exactly.

- ...:

  Absorbed and ignored - lets this function sit behind
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)'s
  `...` passthrough (which also forwards `n_samples`/`grid_resolution`
  meant for other routes) without erroring on arguments this route
  doesn't use.
