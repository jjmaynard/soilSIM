# Apply GP Depth Trends with Correlation Preservation

Enhanced core function that applies GP-derived depth trends using Module
8 utilities for robust error handling and validation.

## Usage

``` r
apply_gp_depth_trends(
  cokey_data,
  gp_predictions,
  properties,
  preserve_correlations = TRUE,
  primary_property = NULL,
  verbose = getOption("ssurgo.verbose", FALSE),
  config = NULL,
  gp_models = NULL
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- gp_predictions:

  Named list of GP predictions by property

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- primary_property:

  Reference property for correlation preservation

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

- config:

  Optional Monte Carlo config (as from
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md))
  whose `monte_carlo$vertical_correlation_method` selects between
  `"joint_copula"` (default as of
  `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md` Phase 13 - dispatches to
  [`preserve_correlation_structure_joint()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure_joint.md),
  drawing depth correlation and property correlation simultaneously) and
  `"gp_quantile_retrofit"` (the original algorithm, still fully
  supported as an explicit opt-out - dispatches to
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)).
  `NULL` (default) resolves to `"joint_copula"`, matching
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md)'s
  own default - set
  `config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit"`
  explicitly to opt back into the original behavior. Under
  `"joint_copula"`, `config$monte_carlo$vertical_correlation_gating`
  (default `FALSE`) separately controls whether `bound_sd`-based
  discontinuity gating (Phase 1c/1d) is applied - kept independent of
  the core method choice since its numeric defaults are not yet
  empirically calibrated (Phase 8).

- gp_models:

  Optional named list of fitted GP models (as
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)
  returns), keyed by property - passed through to
  [`preserve_correlation_structure_joint()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure_joint.md)
  when `vertical_correlation_method = "joint_copula"`, so its depth
  kernel can reuse each property's already-fitted,
  already-cross-validated length-scale (see
  [`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md))
  instead of falling back to a full-depth-range default. Ignored
  entirely under `"gp_quantile_retrofit"`.

## Value

Adjusted simulation data
