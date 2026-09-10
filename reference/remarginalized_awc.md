# Per-Pixel Available Water Capacity from a Re-Marginalized Ensemble

Plant-available water capacity (cm, summed over the ensemble's depth
windows), per realization, then summarized per pixel. Deliberately does
NOT use
[`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)
(ROSETTA + van Genuchten): that needs a live network POST and its own
Monte Carlo per call, so it cannot run per (pixel, realization).

## Usage

``` r
remarginalized_awc(
  mukey_ensemble,
  posterior_by_property_window,
  method = c("saxton_rawls", "direct"),
  n_out = 250,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  rock_fragment = TRUE,
  soc_to_om = 1.724,
  tile_rows = NULL,
  restriction_depth = NULL
)
```

## Arguments

- mukey_ensemble:

  An
  [`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
  result.

- posterior_by_property_window:

  As in
  [`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
  keyed by the ensemble's property (sim-column) names. `"saxton_rawls"`
  needs `sand_total`/`silt_total`/`clay_total`/`db` (+ optional `soc`,
  `rfv`); `"direct"` needs `wr_3b`/`wr_15b` (+ optional `rfv`).

- method:

  `"saxton_rawls"` (default) or `"direct"` - see Method.

- n_out, probs:

  As in
  [`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md).

- rock_fragment:

  If `TRUE` (default) apply the rock-fragment correction when an `rfv`
  posterior is available.

- soc_to_om:

  Multiplier converting the `soc` posterior (organic carbon %) to
  organic matter % for Saxton-Rawls (default 1.724, the Van Bemmelen
  factor). Only used by `"saxton_rawls"`.

- tile_rows:

  Optional. Process the AOI in horizontal row-strips, releasing each
  strip's realization stack before the next. Per-pixel quantiles are
  cell-independent, so a tiled run is numerically identical to
  `tile_rows = NULL`. `NULL` (default) / `>= nrow` disables tiling.

- restriction_depth:

  Optional single-layer
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of restriction depth in cm - typically
  `fetch_solus_restriction_depth(aoi_vect)`, already
  right-censoring-guarded (cells with no detected restriction are `Inf`,
  never `SOLUS_RESTRICTION_CENSOR_CM` itself). `NULL` (default)
  preserves today's behavior exactly - every window contributes its full
  nominal thickness, regardless of bedrock. When supplied, resampled
  once (bilinear - a continuous depth value, not a category) onto the
  working grid, and each window's contribution is scaled by its
  **effective** thickness
  `pmax(0, pmin(bottom, restriction_depth) - top)` instead of the
  nominal `bottom - top`: a window entirely above the restriction keeps
  its full thickness, one straddling it gets partial credit, one
  entirely below it contributes exactly **0** (not `NA` - `NA` would
  incorrectly blank the whole pixel's AWC, including valid shallower
  windows, rather than just this window's contribution). A pixel with
  `restriction_depth` missing/`NA` (e.g. no coverage at the AOI edge) is
  treated the same as `Inf` - no truncation - the same "no evidence -\>
  don't truncate" convention this pipeline already uses for other
  SSURGO/SOLUS coverage gaps, rather than propagating `NA` into the AWC
  total.

## Value

`list(awc_cm = <named list of P.. SpatRasters>, n_kept =, windows_used =, n_tiles =)`.
AWC clamped at 0.

## Method

- **`"saxton_rawls"`** (default):
  `SOLUS100 publishes no water-retention variable`, so `wr_3b`/ `wr_15b`
  can never be fused directly. Instead the fusable
  `sand_total`/`silt_total`/`clay_total`/`db` (plus optional `soc`,
  `rfv`) are re-marginalized to their per-pixel fused posteriors, then
  [`saxton_rawls_raster()`](https://jjmaynard.github.io/soilSIM/reference/saxton_rawls_raster.md)
  derives field capacity / wilting point per (pixel, realization).
  `AWC = sum_windows (fc - wp)/100 * thickness` (the Saxton-Rawls
  conversion already applies the rock-fragment correction internally).

- **`"direct"`**: uses `wr_3b`/`wr_15b` posteriors supplied in
  `posterior_by_property_window` directly
  (`AWC = sum_windows (wr_3b - wr_15b)/100 * thickness * (1 - rfv/100)`).
  Only usable if the caller has a water-retention likelihood from some
  non-SOLUS source.

## Interpreting the per-pixel AWC distribution

`AWC` is computed **per source realization** and combined by
[`terra::quantile()`](https://rspatial.github.io/terra/reference/quantile.html)
per pixel, so `awc_cm$P95 - awc_cm$P5` is a genuine 90% credible
interval. Because
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md)
threads one realization index through every property and every window
(rank-preserving), each realization's inputs are cross-property *and*
vertically consistent - a "wet" profile is wet at all depths. What the
band does and does not capture:

- **Marginal** uncertainty of each input property: the full SSURGO x
  SOLUS fused posterior spread per pixel, propagated through
  Saxton-Rawls. Captured.

- **Dependence** (cross-property, cross-depth): the SSURGO tabular
  copula (KSSL + joint-copula vertical correlation), *held fixed*. SOLUS
  carries no joint information, so nothing is re-estimated. The
  preserved correlation is **rank (Spearman)**, not Pearson.

- **Pedotransfer-function error**: NOT propagated. Saxton-Rawls is
  applied deterministically per realization; the +/-15\\

- **Sub-mukey variation in dependence**: none. Every pixel in a mukey
  shares the source ensemble; only the fused marginals vary
  pixel-to-pixel.

See `docs/09_multi_source_raster_fusion_pipeline.md` -\> "Statistical
structure of the per-pixel bridge" for the full treatment.

## See also

[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`saxton_rawls_raster()`](https://jjmaynard.github.io/soilSIM/reference/saxton_rawls_raster.md),
[`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md),
[`fetch_solus_restriction_depth()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_restriction_depth.md)
