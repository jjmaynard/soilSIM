# Vectorized per-chunk core of `fuse_texture_group()`

Computes, for every cell (row) in `row_mat` at once, the same
prior/likelihood ILR moments + bivariate-normal fusion that
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md) +
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)
compute one cell at a time - see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
PERF comment at its call site.

## Usage

``` r
fuse_texture_group_batch(
  row_mat,
  clay_id,
  sand_id,
  silt_id,
  prior_z,
  lik_z,
  n_mc = 2000,
  max_cells_per_subchunk = 2000
)

fuse_texture_group_batch_core(
  row_mat,
  clay_id,
  sand_id,
  silt_id,
  prior_z,
  lik_z,
  n_mc = 2000
)
```

## Arguments

- row_mat:

  A matrix with one row per cell and the 18
  `<id>_<prior|lik>_<lo|p50|hi>` columns
  [`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
  builds.

- clay_id, sand_id, silt_id:

  The three members' `id`s, in ILR order.

- prior_z, lik_z:

  Standard-normal quantiles for the prior/likelihood low-high intervals
  (as in
  [`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)'s
  `z` argument).

- n_mc:

  Monte Carlo sample size per side (matches
  [`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)'s
  default).

- max_cells_per_subchunk:

  Upper bound on cells processed by one
  `fuse_texture_group_batch_core()` call - each cell needs `6 * n_mc`
  doubles of intermediate MC-draw storage (`n_mc=2000` -\> ~96KB/cell),
  so
  [`terra::app()`](https://rspatial.github.io/terra/reference/app.html)'s
  own chunk size (which only accounts for the raw, non-MC-expanded
  raster layers) can hand
  [`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
  chunks large enough to exhaust memory - confirmed empirically: a
  50,000-cell chunk failed with "cannot allocate vector of size 4.5 Gb".
  Sub-chunking here, in cell order, doesn't change the RNG stream order
  (still cell-major) so results are unaffected.

## Value

A `nrow(row_mat) x 5` matrix, columns `mu1, mu2, S11, S12, S22`.

## Details

RNG draws are generated as a single
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) call ordered so the
underlying standard-normal stream is consumed in exactly the same
cell-major, then prior-before-lik, then clay/sand/silt, then
draw-1..n_mc order the original per-cell loop consumed it in (rnorm()'s
`mean`/`sd` args only affine-transform each already-drawn standard
normal; they don't change how many values are drawn or in what order) -
so this is a pure vectorization, not a behavior change.
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)'s
2x2 [`solve()`](https://rdrr.io/r/base/solve.html)-based fusion is
replaced with the closed-form 2x2 matrix inverse formula
([`solve()`](https://rdrr.io/r/base/solve.html) isn't vectorizable
across cells), which can differ from
[`solve()`](https://rdrr.io/r/base/solve.html) by floating-point
rounding (~1e-10 relative) but is otherwise the identical linear
algebra.
