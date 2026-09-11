# Per-Pixel Ensemble Re-Marginalization (raster-fusion \<-\> tabular MC bridge)

Transforms the per-mukey joint Monte Carlo ensemble from
[`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
(`R/adapter-ssurgo-simulate.R`) onto the SOLUS grid so that each pixel's
marginal distributions match
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)'s
SOLUS-fused posterior, while preserving the ensemble's empirical copula
(cross-property and cross-depth-window rank structure). This is the
"per-pixel" answer to wiring raster fusion into the modelling framework.

[`zonal_distribution_from_posterior()`](https://jjmaynard.github.io/soilSIM/reference/zonal_distribution_from_posterior.md)
is the cheap mukey-collapsed comparison arm (feeds
`observed_data_by_mukey` in
[`simulate_monte_carlo()`](https://jjmaynard.github.io/soilSIM/reference/simulate_monte_carlo.md));
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md)
is the per-pixel transform;
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md)
is the first derived-quantity consumer (available water capacity, via
[`saxton_rawls_raster()`](https://jjmaynard.github.io/soilSIM/reference/saxton_rawls_raster.md)
since SOLUS100 has no water-retention variable). Nothing here touches
either existing pipeline - purely additive.
