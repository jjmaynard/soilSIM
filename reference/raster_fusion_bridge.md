# Per-Pixel Ensemble Re-Marginalization (raster-fusion \<-\> tabular MC bridge)

Transforms the per-mukey joint Monte Carlo ensemble from
[`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
(`R/ssurgo-simulation.R`) onto the SOLUS grid so that each pixel's
marginal distributions match
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)'s
SOLUS-fused posterior, while preserving the ensemble's empirical copula
(cross-property and cross-depth-window rank structure). This is the
"per-pixel" answer to wiring raster fusion into the modelling framework.

[`zonal_distribution_from_posterior()`](https://jjmaynard.github.io/soilSIM/reference/zonal_distribution_from_posterior.md)
is the cheap mukey-collapsed comparison arm (feeds
`observed_data_by_mukey` in
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md));
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md)
is the per-pixel transform;
[`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)
is the first derived-quantity consumer (available water capacity, via
[`saxton_rawls_raster()`](https://jjmaynard.github.io/soilSIM/reference/saxton_rawls_raster.md)
since SOLUS100 has no water-retention variable). Nothing here touches
either existing pipeline - purely additive.
