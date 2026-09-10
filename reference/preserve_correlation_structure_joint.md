# Preserve Correlation Structure via a Joint Depth x Property Copula

Drop-in alternative to
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md) -
same required parameters, in the same order, so existing call sites work
unchanged - that replaces its sequential
`gp_ratio`/single-"primary-property" retrofit with the joint
Kronecker-copula sampler from the joint Kronecker-copula sampler
(sample_joint_depth_property_copula() +
[`apply_copula_to_marginals()`](https://jjmaynard.github.io/soilSIM/reference/apply_copula_to_marginals.md)),
so depth correlation and property correlation are satisfied
SIMULTANEOUSLY by construction rather than approximated by a
rank-copying retrofit. `primary_property` is accepted (for signature
compatibility with
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md))
but unused - the joint method has no need for a single reference
property, since every property's correlation to every other is modeled
directly via `R_prop`.

## Usage

``` r
preserve_correlation_structure_joint(
  property_matrices,
  gp_predictions,
  depths,
  primary_property,
  verbose = getOption("ssurgo.verbose", FALSE),
  gp_models = NULL,
  boundary_distinctness = NULL,
  kernel = c("exponential", "matern")
)
```

## Arguments

- property_matrices:

  Named list of property matrices (rows = depths, columns =
  simulations) - same shape
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)
  consumes.

- gp_predictions:

  Named list of GP-predicted per-depth mean vectors - passed through to
  [`apply_copula_to_marginals()`](https://jjmaynard.github.io/soilSIM/reference/apply_copula_to_marginals.md)'s
  optional GP-mean recentering.

- depths:

  Depth vector (real units, e.g. cm).

- primary_property:

  Accepted for signature compatibility with
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md);
  unused by the joint method.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

- gp_models:

  Optional named list of fitted GP models (as returned by
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md),
  or raw `GPfit`-classed objects), keyed by property. When supplied, the
  depth kernel's length-scale is derived from
  [`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md)
  applied to every property with a usable model, averaged across them
  (reusing each property's already-fitted, already-cross-validated GP
  fit). When `NULL`/empty (e.g. this function is called standalone,
  without the fitted models available), falls back to a length-scale
  spanning the full depth range (`diff(range(depths))`) - a conservative
  "moderate smooth correlation across the whole profile" default that
  keeps this function usable without requiring GP models.

- boundary_distinctness:

  Optional per-depth `bound_sd` vector, passed through to
  [`build_depth_correlation_kernel()`](https://jjmaynard.github.io/soilSIM/reference/build_depth_correlation_kernel.md)
  for discontinuity gating. `NULL` (default) skips gating.

- kernel:

  `"exponential"` (default) or `"matern"` - passed through to
  [`build_depth_correlation_kernel()`](https://jjmaynard.github.io/soilSIM/reference/build_depth_correlation_kernel.md).

## Value

List of adjusted property matrices - same shape/contract as
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)'s
return value. Degrades gracefully to the unadjusted `property_matrices`
(with a warning) on insufficient dimensions or a sampling/ mapping
failure, matching
[`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)'s
own graceful-failure contract.
