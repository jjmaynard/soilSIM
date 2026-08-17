# Component-Composition and Correlated-Triangular Property Simulation

Component-composition simulation
([`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)),
correlated triangular-distribution sampling
([`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)),
and per-cokey flexible-property simulation
([`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)),
ported from `code_ref/brdf/property_simulation.R`. Also includes two
small standalone horizon-data helpers from the same source file
([`remove_organic_layer()`](https://jjmaynard.github.io/soilSIM/reference/remove_organic_layer.md),
[`slice_and_aggregate_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/slice_and_aggregate_soil_data.md)).

Reimplemented to avoid two new dependencies, consistent with this
project's established preference for a small reimplementation over a new
package (already avoided `rmetalog`, `compositions` elsewhere):

- [`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)'s
  uncorrelated-normal draw used
  `MASS::mvrnorm(n, mu = rep(0, k), Sigma = diag(k))` - mathematically
  identical to `k` independent `stats::rnorm(n)` draws when `Sigma` is
  the identity matrix, so no `MASS` dependency is needed.

- [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
  texture step used `compositions::acomp()`/`ilr()`/ `ilrInv()` -
  replaced with `R/distributions.R`'s already-validated
  [`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)/
  [`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)
  (documented there as matching `compositions::ilr()` exactly), exactly
  as the raster-fusion port (`R/raster-fusion.R`) already did for the
  same reason.

`simulate_cokey` (the file's own header flags it as "an earlier,
superseded implementation... may be dead code") and
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
(a name collision with the unrelated, already-real
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
in `R/gp-modeling.R`; duplicates work `R/monte-carlo.R`'s
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)
already does with a better, percentile-triplet-family fitting engine
rather than triangular-only sampling) are deliberately not ported.
