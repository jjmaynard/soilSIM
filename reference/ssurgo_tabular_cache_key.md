# Build a cache key for one AOI's raw SSURGO tabular download (depth-independent)

[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
(`R/ssurgo-acquisition.R`) takes no depth-window argument at all - it
fetches every horizon for every AOI mukey; depth filtering happens
later, per call, in
[`aggregate_depth_window_by_replicate()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_depth_window_by_replicate.md).
The downloaded tabular data therefore depends only on the AOI (and the
fixed default `properties` list every call site uses - never varied),
not on any depth window. Keying the tabular cache on
`top_depth`/`bottom_depth` (as it was before this function existed -
MULTI_PROPERTY_FUSION_PLAN.md task L1) fragments one AOI's cache into
one unusable copy per distinct window ever requested for it, so e.g.
[`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)'s
wide-span simulation and a later single-window
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
call for the same AOI never shared a tabular download. Mirrors
[`mukey_grid_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/mukey_grid_cache_key.md)'s
fixed-sentinel pattern.

## Usage

``` r
ssurgo_tabular_cache_key(aoi_vect)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (single-feature AOI).

## Value

A character string, safe for use as a filename.
