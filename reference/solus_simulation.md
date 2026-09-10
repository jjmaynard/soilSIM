# Raster SOLUS Percentile Likelihood

The SOLUS100 half of the raster fusion prior/likelihood pipeline (see
`R/raster-fusion.R`): fetches SOLUS100 low/prediction/high rasters via
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
for a requested depth window, in the `list(values=, probs=)` shape
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
expects.

[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)'s
output-layer-naming assumption
(`paste0(solus_variable, "_", depth_slice, "_cm_", suffix)`) is verified
against a live
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
call for `output_type` `"prediction"` / `"95% low prediction interval"`
/ `"95% high prediction interval"` (suffixes `p`/`l`/`h`).

`solus_variable` is a direct per-call argument (there is no global
property registry). Note that `awc` is not a valid `fetchSOLUS()`
variable.
