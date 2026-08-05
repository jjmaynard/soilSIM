# Vectorized feasibility check: a valid metalog quantile function must be monotonically increasing in y (equivalent to its density staying non-negative everywhere). Probes `quantile_metalog_linear_raster()` at a fixed grid of y-values and flags cells where consecutive probe values decrease - a probe, not a proof, but fully vectorized raster arithmetic. Streams one probe raster at a time rather than materializing all of them (validated upstream to avoid an allocation failure at large cell counts).

Vectorized feasibility check: a valid metalog quantile function must be
monotonically increasing in y (equivalent to its density staying
non-negative everywhere). Probes
[`quantile_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
at a fixed grid of y-values and flags cells where consecutive probe
values decrease - a probe, not a proof, but fully vectorized raster
arithmetic. Streams one probe raster at a time rather than materializing
all of them (validated upstream to avoid an allocation failure at large
cell counts).

## Usage

``` r
check_metalog_feasibility_raster(
  fit,
  bounds,
  boundedness,
  y_grid = seq(0.02, 0.98, by = 0.02)
)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md).

- bounds, boundedness:

  Same as
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)'s
  arguments for this `fit`.

- y_grid:

  Probability grid to probe.
