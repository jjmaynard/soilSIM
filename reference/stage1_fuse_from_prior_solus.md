# Fuse an Already-Fetched Prior and Likelihood (Stage 1 fusion tail)

The portion of
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
that runs once its SSURGO prior and SOLUS likelihood percentile lists
are in hand: resample the prior onto the SOLUS grid, optionally attach
the real per-mukey Monte Carlo draws for raw-draws fusion, fuse via
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md),
and assemble the standard
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
return list. Factored out so
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
and
[`run_fusion_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_multiproperty.md)
share exactly one non-compositional fusion code path.

## Usage

``` r
stage1_fuse_from_prior_solus(
  property_config,
  prior,
  solus,
  draws = NULL,
  mukey_raster_native = NULL,
  verbose = TRUE
)
```

## Arguments

- property_config:

  As in
  [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md).

- prior, solus:

  The pre-alignment `list(values=, probs=)` percentile lists for the
  SSURGO prior and SOLUS likelihood.

- draws:

  Optional in-memory draws data frame (already simulated by the caller)
  for the raw-draws fusion path; `NULL` skips it.

- mukey_raster_native:

  Optional native-grid mukey raster matching `draws` - required for the
  raw-draws path.

- verbose:

  Passed through to
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
  (default `TRUE`, matching the original inline behavior).

## Value

[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)'s
return list
(`list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=, route_detail=, n_fallback_cells=)`).
