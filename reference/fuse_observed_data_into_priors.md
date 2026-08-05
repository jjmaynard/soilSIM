# Fuse Observed Field/Lab Data into SSURGO-Derived Priors

Orchestrates Bayesian updating against each property's per-horizon prior
fit produced by
[`prepare_simulation_parameters()`](https://jjmaynard.github.io/soilSIM/reference/prepare_simulation_parameters.md),
using `bayesian-updating.R`'s pure fusion primitives
([`bayes_fuse()`](https://jjmaynard.github.io/soilSIM/reference/bayes_fuse.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)).
For every property present in `observed_data`, every horizon's existing
prior is replaced by its fused posterior - each horizon keeps its OWN
prior (so posteriors still vary per horizon), but all incorporate the
SAME shared observed-data likelihood, since field/lab measurements
aren't typically tied to individual simulated horizons at that
granularity.

## Usage

``` r
fuse_observed_data_into_priors(
  simulation_params,
  simulation_data,
  sim_properties,
  properties,
  observed_data,
  composition_plan,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_params:

  Per-horizon parameter list from
  [`prepare_simulation_parameters()`](https://jjmaynard.github.io/soilSIM/reference/prepare_simulation_parameters.md).

- simulation_data:

  The filtered per-horizon data frame (same row order as
  `simulation_params`).

- sim_properties:

  Character vector `simulation_params` is keyed over (may include
  composition-group pseudo-properties).

- properties:

  The caller-facing, original property vector. Unused directly (fusion
  is keyed by `sim_properties`/`observed_data` names); kept for
  interface symmetry with sibling pipeline functions.

- observed_data:

  Named list keyed by property name (see shapes above).

- composition_plan:

  Result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md).

- config:

  Simulation configuration.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

`simulation_params`, with fused entries replaced in place.

## Details

Two supported shapes per non-texture entry in `observed_data`:

- a numeric vector of raw samples -\> the fully general route
  ([`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md)'s
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
  grid-KDE path). The prior is first SAMPLED (via
  [`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md),
  since that route needs both sides to already be raw vectors - it fits
  no parametric family to either side), and the resulting posterior
  sample is stored as an exact empirical quantile function
  (`family = "linear_cdf"`) rather than forced back into a possibly
  ill-fitting parametric family.

- a family-native parameter list (`list(mean=, sd=)` for a normal or
  lognormal prior/likelihood; `list(shape1=, shape2=)` or
  `list(mean=, var=)` for a beta prior/likelihood) -\> the closed-form
  route
  ([`bayes_fuse()`](https://jjmaynard.github.io/soilSIM/reference/bayes_fuse.md)),
  only supported for prior families `"normal"`/`"lognormal"`/`"beta"`
  (the only families with a conjugate route). Any other prior family
  (`"triangular"`/`"uniform"`/ `"metalog"`/`"linear_cdf"`) SKIPS fusion
  for that (horizon, property) with a logged WARN, continuing to
  simulate from the unfused prior rather than erroring the whole
  pipeline.

The texture composition group (if active) is handled jointly: supply all
three of `claytotal`/`sandtotal`/`silttotal` in `observed_data` (each
either a `c(low=,rep=,high=)` triplet, an unnamed length-3
`c(low,rep,high)` vector, or a raw sample vector reduced to one via its
own empirical quantiles at `config$monte_carlo$lh_percentile`) - partial
texture entries (1-2 of 3) skip texture fusion entirely with a WARN,
matching the same "some but not all composition members present"
fallback
[`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md)
already uses. Fusion runs once per horizon via
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)
(that horizon's own clay/sand/silt l/r/h as the prior side, the shared
observed triplet as the likelihood side), replacing
`simulation_params[[i]][["ilr1"]]`/`[["ilr2"]]`.

Note: this only narrows each property's MARGINAL prior spread. It does
NOT change the cross-property correlation/dependency structure used by
the Cholesky copula in
[`simulate_correlated_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_properties.md) -
[`estimate_property_correlations()`](https://jjmaynard.github.io/soilSIM/reference/estimate_property_correlations.md)
(Step 5, when `auto_correlation=TRUE`) estimates that from
`simulation_data`'s raw representative values, which fusion does not
touch.
