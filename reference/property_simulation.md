# Component-Composition and Correlated-Triangular Property Simulation

Component-composition simulation
([`simulate_component_composition()`](https://jjmaynard.github.io/soilSIM/reference/simulate_component_composition.md)),
correlated triangular-distribution sampling
([`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)),
and per-cokey flexible-property simulation
([`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)),
plus two small standalone horizon-data helpers
([`remove_organic_layer()`](https://jjmaynard.github.io/soilSIM/reference/remove_organic_layer.md),
[`slice_and_aggregate_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/slice_and_aggregate_soil_data.md)).

[`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)
draws uncorrelated normals with
[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html) (no `MASS`
dependency), and
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
texture step uses `R/core-distributions.R`'s
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)
(no `compositions` dependency).
