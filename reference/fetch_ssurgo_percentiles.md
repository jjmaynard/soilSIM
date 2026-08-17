# Fetch SSURGO Percentile-Value Rasters for an AOI

The top-level SSURGO "prior" entry point for `R/raster-fusion.R`'s
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md):
rasterizes map units, Monte Carlo-simulates the requested property, and
returns per-cell percentile-value rasters in the `list(values=, probs=)`
shape
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
expects.

## Usage

``` r
fetch_ssurgo_percentiles(
  aoi_vect,
  property_id,
  top_depth,
  bottom_depth,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  n_mc = 1000,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- property_id:

  One of
  [`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md)'s
  recognized ids.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- probs:

  Percentile probabilities to compute (default
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`).

- n_mc:

  Number of triangular draws per component (default 1000).

- parallel, n_cores:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
  `parallel`/`n_cores` - the per-cokey depth-trend GP fitting step is
  this function's dominant cost for AOIs with many cokeys. Default
  `parallel = FALSE` matches prior behavior exactly.

## Value

`list(values = <named list of percentile-value SpatRasters>, probs = probs)`,
or `NULL` if the mukey raster or the Monte Carlo draws are unavailable
for this AOI.
