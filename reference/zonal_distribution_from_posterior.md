# Zonal (per-mukey) Reduction of a Fused Posterior Raster - Benchmark Baseline

Collapses
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)'s
per-pixel posterior to one distribution per mukey by taking
[`terra::zonal()`](https://rspatial.github.io/terra/reference/zonal.html)
means of its percentile layers, in the `low`/`rep`/`high` shape
[`fuse_observed_data_into_priors()`](https://jjmaynard.github.io/soilSIM/reference/fuse_observed_data_into_priors.md)'s
`observed_data_by_mukey` accepts. **Percentile-based only - there is
deliberately no point-estimate ("mean") path**, which would narrow the
fed-back prior and understate uncertainty (see the design doc's Problem
A redraft).

## Usage

``` r
zonal_distribution_from_posterior(
  posterior,
  mukey_raster,
  probs = c(0.05, 0.5, 0.95)
)
```

## Arguments

- posterior:

  A
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)-style
  `list(percentiles = , ...)`.

- mukey_raster:

  A categorical mukey
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  (any grid - resampled to the posterior grid via nearest-neighbour
  internally).

- probs:

  Length-3 ascending probabilities for the `low`/`rep`/`high` triplet;
  each must be one of the posterior's own percentile-layer probabilities
  (no interpolation).

## Value

Named list keyed by mukey code, each element `c(low =, rep =, high =)`.

## Details

Its sole purpose is the A.6 benchmark: the mukey-collapsed comparison
arm for
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
run through the *same* tabular Monte Carlo machinery via
`observed_data_by_mukey`. Kept `@keywords internal` until that benchmark
says whether the cheap path is adequate for production (then promote) or
not (then delete).
