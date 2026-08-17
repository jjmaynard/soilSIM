# Fuse one property's prior and likelihood, dispatching to the right family/route automatically from `property_config`.

The raster-native counterpart of `bayesian-updating.R`'s
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
and this toolkit's top-level entry point for non-compositional
properties (see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
for the compositional/texture case) - the raster analogue of
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
minus the SSURGO/SOLUS fetch-and-cache wrapper (see this file's
`@section Deliberately out of scope:` above).

## Usage

``` r
fuse_property_adaptive(
  prior_value_rasters,
  prior_probs,
  lik_value_rasters,
  lik_probs,
  property_config,
  threshold_cells = 80000,
  ...
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Named lists of percentile-value SpatRasters.

- prior_probs, lik_probs:

  Matching probability vectors (reconciled via
  [`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md)
  if they differ).

- property_config:

  List with `dist` (one of "normal"/"beta"/"gamma"/
  "lognormal"/"metalog"/"auto"), and optionally `bounds`/`boundedness`/
  `auto_skew_threshold`.

- threshold_cells:

  AOI cell count at or below which
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)
  uses its general route.

- ...:

  Additional arguments passed to
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)/
  `fuse_lognormal_adaptive()` (e.g. `n_samples`, `grid_resolution`,
  `verbose`).

## Value

`list(posterior=, route=, route_detail=, n_fallback_cells=, dist=, dist_source=, skew_proxy=)`.
