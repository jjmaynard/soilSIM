# Fetch SOLUS100 Low/Prediction/High Rasters for Multiple Variables, One Depth Window (S1)

Batched sibling of
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md):
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
accepts `variables`, `depth_slices`, and `output_type` all as vectors
simultaneously (live-verified 2026-09-04 - see ), so every variable's
low/prediction/high rasters for a window can be fetched in **one**
network call instead of `3 * length(solus_variables)` separate ones. The
returned layer-naming convention
(`paste0(variable, "_", depth, "_cm_", suffix)`) was confirmed identical
to the single-variable case for the fully-batched request (`nlyr`
matched `length(variables) * length(depth_slices) * length(output_type)`
exactly, no layers dropped,
[`names()`](https://rdrr.io/r/base/names.html) order not request order
but that does not matter since layers are looked up by exact name).

## Usage

``` r
fetch_solus_low_pred_high_multiproperty(
  aoi_vect,
  solus_variables,
  top_depth,
  bottom_depth
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- solus_variables:

  Character vector of
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable names.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

## Value

A named list, one entry per element of `solus_variables`, each
`list(pred=, low=, high=)` in the same shape
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
returns (a single-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
or `NULL` per output type). If the single underlying `fetchSOLUS()` call
itself errors, every variable's entry is
`list(pred=NULL, low=NULL, high=NULL)` (a whole-request failure, not a
crash) - callers should fall back to the scalar
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
per variable in that case. A per-`(variable, output_type)` combo missing
from an otherwise-successful response is `NULL` for just that piece,
with a [`warning()`](https://rdrr.io/r/base/warning.html) naming it - a
partial failure, not a whole- request one.

## See also

[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md),
[`fetch_solus_percentiles_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles_multiproperty.md)
