# Raster SOLUS Percentile Likelihood

The SOLUS100 half of the raster fusion prior/likelihood pipeline (see
`R/raster-fusion.R`): fetches SOLUS100 low/prediction/high rasters via
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
for a requested depth window, in the `list(values=, probs=)` shape
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
expects. Ported from `code_ref/reanalysis-platform/solus_raster.R`.

[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)'s
output-layer-naming assumption
(`paste0(solus_variable, "_", depth_slice, "_cm_", suffix)`) was
verified against a live, network-connected
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
call before porting (not assumed from the source comments alone) -
confirmed exact for `output_type`
`"prediction"`/`"95% low prediction interval"`/`"95% high prediction interval"`
(suffixes `p`/`l`/`h`).

Unlike the original bundle's `config.R`-driven `solus_variable` lookup
(which had a confirmed-wrong `awc = "awc"` entry - not a valid
`fetchSOLUS()` variable), `solus_variable` is a direct per-call argument
here, consistent with `R/raster-fusion.R`'s decision not to port the
`PROPERTIES`/`config.R` global registry.
