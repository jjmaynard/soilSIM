# Fuse one property's prior and likelihood, dispatching to the right family/route automatically from `property_config`.

The raster-native counterpart of `core-fusion.R`'s
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
and this toolkit's top-level entry point for non-compositional
properties (see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
for the compositional/texture case) - the raster analogue of
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
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
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)/[`fuse_lognormal_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_lognormal_adaptive.md)/
  [`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)
  (e.g. `n_samples`, `grid_resolution`, `verbose`, `posterior_probs` -
  see
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md);
  reaches every route -
  [`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)
  absorbs and ignores the arguments it doesn't use, e.g. `n_samples`).
  Also how `mukey_raster`/`mukey_draws` (opting into
  `prior_fusion_method = "raw_draws"`, see
  [`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
  docs) reach the underlying routes: both the general-KDE and
  closed-form routes support it for
  `dist %in% c("normal","beta","gamma")`, and it reaches
  [`fuse_lognormal_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_lognormal_adaptive.md)'s
  small-AOI general branch AND large-AOI closed-form branch the same
  way.
  [`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)
  also has a raw-draws fit (empirical per-mukey percentiles feed the
  same interpolation machinery) - but
  [`resolve_want_raw_draws()`](https://jjmaynard.github.io/soilSIM/reference/resolve_want_raw_draws.md)'s
  default still excludes `dist = "metalog"`, mirroring the texture-group
  route's treatment: a fit existing is not the same as its cost/benefit
  being benchmarked, and that default decision is kept separate.
  Explicitly passing `mukey_raster`/`mukey_draws` for `dist = "metalog"`
  still works.

## Value

`list(posterior=, route=, route_detail=, n_fallback_cells=, dist=, dist_source=, skew_proxy=)`.
