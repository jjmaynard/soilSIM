# Raster SSURGO Percentile Prior: Mukey Lookup and Monte Carlo Draws

The SSURGO half of the raster fusion prior/likelihood pipeline (see
`R/raster-fusion.R`): given a
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
AOI, rasterizes SSURGO map units, Monte Carlo-simulates per-component
soil properties (reusing already-ported
[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)/[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
from `R/property-simulation.R`), aggregates to a requested depth window,
and rasterizes per-mukey percentiles for
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
to consume. Ported from
`code_ref/reanalysis-platform/{ssurgo_simulation.R, ssurgo_prior.R}`.

Reuses `R/kssl-reference-correlations.R`'s already-built-in,
already-shipped
[`.kssl_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/dot-kssl_property_matrices.md)/[`.kssl_texture_matrices()`](https://jjmaynard.github.io/soilSIM/reference/dot-kssl_texture_matrices.md)
genhz-keyed correlation matrices directly as
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
`correlation_matrices`/ `txt_correlation_matrices` arguments, rather
than re-reading the raw
`data/global_cor_matrices.rds`/`data/global_cor_texture_matrices.rds`
files at runtime the way the source's
`load_global_correlation_matrices()` did (those files already back
`R/sysdata.rda`, built once via
`data-raw/build_kssl_reference_correlations.R` - no new asset-loading
code needed). Similarly reuses
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
(already exported from `kssl-reference-correlations.R`) in place of the
source's
[`aqp::generalizeHz()`](https://ncss-tech.github.io/aqp/reference/generalize.hz.html)
call.

[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
(`R/ssurgo-acquisition.R`) default `properties` and always-included base
columns (`comppct_l/r/h`, all `_l/_r/_h` horizon triplets) already match
[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)/[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
expected SSURGO-stem input vocabulary directly - the only translation
table genuinely needed is
[`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md),
mapping a caller-facing property id to
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
short-code *output* column name.
