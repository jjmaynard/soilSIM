# Align a Prior and Likelihood onto One Shared Grid

Every fusion route requires cell-aligned prior/likelihood rasters, so
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md),
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md),
and
[`run_fusion_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_multiproperty.md)
each resample one side onto the other's grid before fusing. Factored out
here so the resampling *direction* is one source-agnostic decision point
(prior vs. likelihood, not "SSURGO" vs. "SOLUS") rather than several
copies of the same
[`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html)
call.

## Usage

``` r
align_prior_likelihood(prior_values, lik_values, resampling = c("down", "up"))
```

## Arguments

- prior_values, lik_values:

  Named lists of percentile `SpatRaster`s (a `prior$values` /
  `likelihood$values`-shaped list, as
  [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
  etc. hold them).

- resampling:

  `"down"` (default, matches all prior behavior exactly) resamples the
  prior DOWN onto the likelihood's grid via bilinear interpolation - the
  fused output sits at the likelihood's own resolution, which for the
  SSURGO/SOLUS100 pipeline means SSURGO's finer native grid (~30 m) is
  blurred/averaged into SOLUS100's coarser 100 m cells.

  `"up"` instead resamples the likelihood UP onto the prior's grid via
  nearest-neighbor ("block replication" - every fine cell takes its
  enclosing coarse cell's exact value, with no invented intermediate
  values) - the fused output sits at the prior's own (finer) resolution
  instead, preserving the prior's real spatial detail (e.g. SSURGO
  map-unit boundaries) rather than discarding it. Nearest-neighbor, not
  bilinear, is deliberate here: bilinear would interpolate BETWEEN two
  independent likelihood-cell predictions, implying a smooth underlying
  gradient the likelihood source never actually asserted - it does not
  add real resolution to the likelihood, only fabricates it. Every fine
  cell's likelihood value under `"up"` is traceable back to one real
  published likelihood cell value.

## Value

`list(prior_values=, lik_values=, reference_grid=)` -
`prior_values`/`lik_values` are both now on `reference_grid` (whichever
input's own grid the other side was resampled onto), returned for reuse
aligning e.g. a categorical mukey raster onto the same grid.
