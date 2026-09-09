# Fetch SOLUS100 Rasters for a Site-Level (Depth-Independent) Variable (S2)

Some
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
variables are **site-level** - a single prediction per pixel with no
depth dependence (e.g. `anylithicdpt`/`resdept`, depth to
bedrock/restriction) - requested via the special `depth_slices = "all"`
value rather than a numeric depth. Confirmed live during the S1 spike
(2026-09-04): the returned layer is named
`paste0(variable, "_all_cm_", suffix)`, exactly one layer per output
type, so unlike
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
there is no depth-window weighting to do - deliberately **not** built on
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md),
which only makes sense for numeric depth slices.

## Usage

``` r
fetch_solus_site_level(
  aoi_vect,
  variable,
  output_types = c("prediction", "95% low prediction interval",
    "95% high prediction interval")
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- variable:

  A
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  site-level variable name (e.g. `"anylithicdpt"`, `"resdept"`).

- output_types:

  Which output types to fetch - any of `"prediction"`,
  `"95% low prediction interval"`, `"95% high prediction interval"`.
  Default fetches all three.

## Value

`list(pred=, low=, high=)`, each a single-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
(or `NULL` for an output type not requested, or not returned by
`fetchSOLUS()`).

## See also

[`fetch_solus_restriction_depth()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_restriction_depth.md)
