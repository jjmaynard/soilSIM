# Extract a Real-Units Depth Length-Scale from a Fitted GP Model

[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) fits a
power-exponential (or Matern) correlation function
`R(h) = corr(0, h; beta, correlation_param)` over depths already
rescaled to `[0,1]` by
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md).
This reuses that fitted correlation function exactly (via
[`GPfit::corr_matrix()`](https://rdrr.io/pkg/GPfit/man/corr_matrix.html),
the same internal machinery `GP_fit()`/`predict.GP()` use) to solve for
the distance `h` (in scaled `[0,1]` units) at which the correlation
first drops to `target_corr` (default `exp(-1)`, the conventional
"range" definition for both exponential and Matern kernels), then
rescales `h` back to real depth units (cm) using the same
`depth_scaling`
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)
stored alongside the model - so no separate/new length-scale estimation
step is introduced; the value is entirely derived from the
already-validated GP fit.

## Usage

``` r
extract_depth_length_scale(
  gp_model_list,
  depth_scaling = NULL,
  target_corr = exp(-1)
)
```

## Arguments

- gp_model_list:

  Either the full list returned by
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)
  (`list(gp_model, depth_scaling, ...)`) or a raw `GPfit`-classed `"GP"`
  object (in which case `depth_scaling` must be supplied separately, or
  the result stays in scaled `[0,1]` units).

- depth_scaling:

  Optional `list(min=, max=, range=)` as stored by
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md);
  only needed when `gp_model_list` is a raw `GP` object without its own
  `$depth_scaling`. If neither is available, the returned length-scale
  stays in scaled `[0,1]` units.

- target_corr:

  The correlation value defining "range" (default `exp(-1) ~= 0.368`).

## Value

A single numeric length-scale in real depth units (or scaled `[0,1]`
units if no `depth_scaling` is available), or `NA_real_` if
`gp_model_list` doesn't contain a usable fitted model (mirrors
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)'s
own `NULL`-on-failure contract - callers should treat `NA_real_` the
same way they'd treat a `NULL` GP model).
