# Default `bayesian_update()` grid resolution used by `fuse_general_kde()` (per-cell KDE fusion route), coarser than `bayesian_update()`'s own standalone default of `0.01`.

Default
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
grid resolution used by
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
(per-cell KDE fusion route), coarser than
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
own standalone default of `0.01`.

## Usage

``` r
FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION
```

## Format

An object of class `numeric` of length 1.

## Why 0.1, not 0.01

[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
[`stats::density()`](https://rdrr.io/r/stats/density.html) calls (kernel
density estimation over `seq(grid_min, grid_max, by = grid_resolution)`)
are the dominant cost of
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
per-cell loop - [`Rprof()`](https://rdrr.io/r/utils/Rprof.html)
profiling on a synthetic 10,000-cell raster
(PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4) attributed 65% of total
wall-clock time to [`density()`](https://rdrr.io/r/stats/density.html)
(`dnorm`/`fft` internals) at the `0.01` default, vs. 12% for the
per-cell percentile-sampling step. At `0.01`, a typical soil-property
range (e.g. 20-50%) produces a ~3,000-point evaluation grid per cell,
per side. Coarsening to `0.1` cuts that grid ~10x and measured **~2.4x**
faster wall-clock on a real
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)
call (23.86s -\> 9.91s at 2,024 cells), while an accuracy sweep across
four representative percentile scenarios (narrow/wide/skewed
distributions) found the posterior mean/variance error at `0.1` stays
under 0.02%/0.12% respectively vs. the `0.01` reference - far inside the
~5-10% variability typical of field-measured soil properties (the same
bar the analogous cubic-spline-to-linear-interpolation tradeoff used
elsewhere in this project's Python sibling report was validated
against). Going coarser (`0.25`+) is where error starts compounding fast
(variance error reaches double digits to 470% by `2.0`) - `0.1` is the
sweet spot, not an arbitrary round number.

[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
own standalone default is deliberately left at `0.01` - this constant
only overrides the grid resolution
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
requests, so any other direct caller of
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
is unaffected.
