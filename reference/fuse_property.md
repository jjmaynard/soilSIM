# Fuse a prior and likelihood belief distribution, dispatching by input shape

Unlike the raster-native reference
([`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)),
which dispatches between the general and closed-form routes by AOI cell
count - a concept with no tabular analogue - this dispatches on the
SHAPE of `prior`/ `likelihood`: atomic numeric vectors (raw samples)
route to the fully general
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md);
named lists of family-native parameters route to the closed-form
[`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md).
This removes the "which heuristic" ambiguity entirely, since the
caller's input shape already determines which route applies.

## Usage

``` r
fuse_property(prior, likelihood, ...)

# Default S3 method
fuse_property(
  prior,
  likelihood,
  family = NULL,
  bounds = NULL,
  method = NULL,
  n_samples = 1000,
  grid_resolution = 0.01,
  ...
)

# S3 method for class 'percentile_rasters'
fuse_property(
  prior,
  likelihood,
  prior_probs,
  lik_probs,
  property_config,
  threshold_cells = 80000,
  ...
)
```

## Arguments

- prior, likelihood:

  Either numeric vectors of raw samples (both sides must be vectors), or
  named lists of family-native parameters matching `family` (both sides
  must be lists).

- ...:

  Passed to the dispatched method.

- family:

  Required when `prior`/`likelihood` are parameter lists; one of
  `"normal"`, `"beta"`, `"gamma"`. Ignored (and inferred as "general")
  when `prior`/`likelihood` are raw-sample vectors.

- bounds:

  Unused currently; kept for interface symmetry with
  percentile-triplet-based callers.

- method:

  Optional assertion/override: `"general"` or `"closed_form"`. Errors if
  it doesn't match what the input shape implies, rather than silently
  overriding it.

- n_samples, grid_resolution:

  Passed to
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
  for the general route.

- prior_probs, lik_probs:

  Percentile levels for `prior`/`likelihood` (raster method only).

- property_config:

  A property configuration list selecting the fusion distribution family
  (raster method only) - see
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

- threshold_cells:

  Raster method only; see
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

## Value

For the general route: a numeric vector of posterior samples
([`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
native output). For the closed-form route: the family-native posterior
parameter list
([`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md)'s
native output). These are NOT the same shape - documented deliberately
rather than forced into a fake-uniform contract, since the caller
already knows which shape it's passing in and therefore which shape it
gets back.

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
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
