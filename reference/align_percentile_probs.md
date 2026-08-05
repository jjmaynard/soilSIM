# Re-express one side's percentile-value rasters onto the other side's probability grid when they differ (e.g. one source's 5th/95th percentiles vs another's 2.5th/97.5th) - every fusion route below fits both sides using ONE shared `probs` vector, so passing mismatched probabilities through unchanged would silently mis-fit whichever side's actual percentiles don't match the assumed labels.

Uses
[`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md)
to reinterpolate the wider-ranging side onto the narrower side's exact
probabilities (the narrower side is always safely interpolatable FROM
the wider one, never the reverse, since
[`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md)
requires the target probability to fall within the source's own range).

## Usage

``` r
align_percentile_probs(
  prior_value_rasters,
  prior_probs,
  lik_value_rasters,
  lik_probs
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Lists of percentile-value SpatRasters.

- prior_probs, lik_probs:

  Matching probability vectors.

## Value

`list(prior_value_rasters=, lik_value_rasters=, probs=)`.
