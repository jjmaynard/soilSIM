# Batched piecewise-linear inverse-CDF sampler across many rows sharing the same `probs` knots

Vectorized equivalent of calling
`sim_linear_cdf(probs, values_mat[i, ], n)` once per row of `values_mat`
and rbind-ing the results - built for
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
(`R/raster-fusion.R`), where
[`stats::approxfun()`](https://rdrr.io/r/stats/approxfun.html) was
previously rebuilt once per raster cell
([`Rprof()`](https://rdrr.io/r/utils/Rprof.html) profiling attributed
~9% of
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
total wall-clock to
[`extract_percentile_pairs()`](https://jjmaynard.github.io/soilSIM/reference/extract_percentile_pairs.md)'s
per-cell dispatch overhead alone, on top of the sampling itself - see
PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4). All rows must share the same
`probs` knots (true for a raster chunk, where every cell's percentile
columns are the same fixed set, e.g. P5/P50/P95) - this is what makes
batching valid; per-row-varying `probs` would need the per-row
[`approxfun()`](https://rdrr.io/r/stats/approxfun.html) approach this
function replaces.

## Usage

``` r
sim_linear_cdf_batch(probs, values_mat, n)
```

## Arguments

- probs:

  Numeric probabilities (0-1), sorted ascending, shared by every row.

- values_mat:

  Numeric matrix, `nrow(values_mat)` rows x `length(probs)` columns
  (column order matching `probs` order), no `NA`/missing values.

- n:

  Number of draws per row.

## Value

A `nrow(values_mat)` x `n` numeric matrix, one row of draws per input
row.

## Details

Matches
[`sim_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf.md)'s
`rule = 2` constant-extrapolation behavior for draws outside
`[min(probs), max(probs)]` (clamped to the first/last `values` column)
exactly, but does NOT preserve
[`sim_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf.md)'s
per-row RNG stream position - the two are statistically equivalent (both
draw `runif(n)` per row and linearly interpolate), not bit-identical,
since this function draws the full `nrow(values_mat) * n` uniform block
in one call rather than one `runif(n)` call per row.
