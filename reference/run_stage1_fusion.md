# Run Stage 1 Fusion for One Property/Depth over an AOI

Fetches the (cached) SSURGO prior (`R/ssurgo-simulation.R`'s
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md))
and SOLUS likelihood (`R/solus-simulation.R`'s
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md))
percentiles, aligns the SSURGO grid onto the SOLUS grid via
[`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html),
then fuses via
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

## Usage

``` r
run_stage1_fusion(
  aoi_vect,
  property_config,
  top_depth,
  bottom_depth,
  composition_groups = NULL,
  property_configs = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI (projected, e.g. EPSG:5070).

- property_config:

  A list with `id` (used as
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)'s
  property id and as the cache key), `solus_variable` (a
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable name), and
  `dist`/`bounds`/`boundedness`/`auto_skew_threshold`/`composition_group`
  as needed by
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

- top_depth, bottom_depth:

  Numeric depth bounds in cm.

- composition_groups, property_configs:

  Only needed when `property_config$composition_group` is set - passed
  through to
  [`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
  (see its docs).

## Value

`list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=, route_detail=, n_fallback_cells=)`,
or `NULL` if the SSURGO or SOLUS side failed. `posterior`'s shape
depends on `dist` - see
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)'s
docs.

## Details

Compositional properties (`property_config$composition_group` set) are
detected up front and delegated to
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
instead, since they must be fetched/fit jointly.
