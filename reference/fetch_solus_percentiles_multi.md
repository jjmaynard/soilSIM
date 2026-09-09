# Fetch SOLUS100 Percentile-Value Rasters for Multiple Variables, One Depth Window (S1)

Batched sibling of
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md):
thin per-variable `list(values=, probs=)` packaging over
[`fetch_solus_low_pred_high_multi()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high_multi.md)'s
one-call fetch.

## Usage

``` r
fetch_solus_percentiles_multi(
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

A named list, one entry per element of `solus_variables`, each either
`list(values = list(P025=, P50=, P975=), probs = c(0.025, 0.5, 0.975))`
(matching
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)'s
shape) or `NULL` if any of that variable's low/pred/high rasters is
unavailable.

## See also

[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md),
[`fetch_solus_low_pred_high_multi()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high_multi.md)
