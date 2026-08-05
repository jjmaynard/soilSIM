# Fuse a prior and likelihood belief distribution, dispatching by input shape

Unlike the raster-native reference
([`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)),
which dispatches between the general and closed-form routes by AOI cell
count - a concept with no tabular analogue - this dispatches on the
SHAPE of `prior`/ `likelihood`: atomic numeric vectors (raw samples)
route to the fully general
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md);
named lists of family-native parameters route to the closed-form
[`bayes_fuse()`](https://jjmaynard.github.io/soilSIM/reference/bayes_fuse.md).
This removes the "which heuristic" ambiguity entirely, since the
caller's input shape already determines which route applies.

## Usage

``` r
fuse_property(
  prior,
  likelihood,
  family = NULL,
  bounds = NULL,
  method = NULL,
  n_samples = 1000,
  grid_resolution = 0.01
)
```

## Arguments

- prior, likelihood:

  Either numeric vectors of raw samples (both sides must be vectors), or
  named lists of family-native parameters matching `family` (both sides
  must be lists).

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
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
  for the general route.

## Value

For the general route: a numeric vector of posterior samples
([`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
native output). For the closed-form route: the family-native posterior
parameter list
([`bayes_fuse()`](https://jjmaynard.github.io/soilSIM/reference/bayes_fuse.md)'s
native output). These are NOT the same shape - documented deliberately
rather than forced into a fake-uniform contract, since the caller
already knows which shape it's passing in and therefore which shape it
gets back.
