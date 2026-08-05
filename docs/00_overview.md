# soilSIM: Architecture Overview

*Entry point for the detailed developer documentation in this directory. Read this first, then
follow the links below into whichever functional area you need.*

## What soilSIM is

`soilSIM` simulates soil properties and fuses estimates from multiple soil data products under
uncertainty. Its core statistical machinery - percentile-triplet and arbitrary-percentile
distribution fitting, correlated Monte Carlo simulation, Gaussian-process depth-trend modeling, and
Bayesian updating/fusion (scalar and raster-native) - is **data-source-agnostic by design**: it
operates on any low/representative/high triplet or prior/likelihood pair, regardless of where the
numbers came from. SSURGO and SOLUS100 are the two concrete data-source **adapters** wired in
today, supplying the prior and likelihood sides (respectively) of the raster fusion pipeline. The
architecture is deliberately split this way so that future data sources (e.g. HWSD, SoilGrids) can
plug in as new adapters alongside SSURGO/SOLUS, reusing the same simulation/fusion core rather than
duplicating it. See `soilSIM/R/soilSIM-package.R` (rendered as `?soilSIM`) for the package-level
roxygen summary of this framing, and [[project-soilsim-package]] in project memory for the decision
history behind it.

This directory's 10 documents go one level deeper than that package doc or the per-function
`.Rd`/roxygen reference: each covers a functional group's full function signatures, internal call
graphs, cross-file dependencies, data flow, and known limitations - the same depth the legacy
`modules/mod0-8_summary.md` docs provided for the pre-package prototype, adapted to soilSIM's
current (larger, restructured) `R/` layout.

## Core vs. adapters

```Shell
                      ┌─────────────────────────────────────────┐
                      │              GENERIC CORE                │
                      │  (data-source-agnostic; works on any     │
                      │   percentile triplet / prior-likelihood  │
                      │   pair, regardless of origin)            │
                      │                                           │
                      │  03 Distribution fitting & correlations  │
                      │  04 Monte Carlo simulation                │
                      │  05 GP depth modeling & adjustment        │
                      │  08 Bayesian updating (scalar)            │
                      │  09 Raster fusion core (distribution-     │
                      │     fitting-raster.R, raster-fusion.R,    │
                      │     raster-cache.R)                       │
                      └───────────────┬───────────────────────────┘
                                      │ consumed by / feeds
              ┌───────────────────────┼───────────────────────┐
              │                       │                       │
   ┌──────────▼──────────┐ ┌──────────▼──────────┐ ┌──────────▼──────────┐
   │  SSURGO ADAPTER      │ │  SOLUS100 ADAPTER    │ │  (future adapters:   │
   │  01 Acquisition/     │ │  09 solus-           │ │   HWSD, SoilGrids,   │
   │     processing        │ │     simulation.R      │ │   etc. - not yet     │
   │  06 Property/depth/   │ │  supplies raster      │ │   implemented)       │
   │     component sim     │ │  fusion likelihood    │ │                      │
   │  09 ssurgo-           │ │  side                 │ │                      │
   │     simulation.R      │ │                       │ │                      │
   │  supplies raster      │ │                       │ │                      │
   │  fusion prior side     │ │                       │ │                      │
   └──────────────────────┘ └──────────────────────┘ └──────────────────────┘
```

`02 Statistics & diagnostics`, `07 AWS/Van Genuchten modeling`, and `10 Utilities` sit alongside
this split: Statistics & Diagnostics validates output from the core and both adapters; AWS modeling
is a self-contained leaf that consumes component/texture data (from either adapter's simulation
output) but nothing calls into it; Utilities is the shared leaf dependency everything else builds
on (validation, logging, config, I/O - see `10_utilities.md` for the one narrow exception where it
reaches back into the SSURGO adapter for a property-lookup helper).

## End-to-end workflow (SSURGO/SOLUS raster fusion path)

```
AOI (WKT) + property list
        │
        ▼
┌───────────────────────────┐
│ 01 download_ssurgo_tabular │──► raw SSURGO tabular data (cached)
└───────────────────────────┘
        │
        ▼
┌───────────────────────────┐
│ 01 process_ssurgo_data /   │──► cleaned, standardized horizon/component tables
│    infill_soil_property    │
└───────────────────────────┘
        │
        ▼
┌───────────────────────────┐
│ 03 fit_percentile_triplet  │──► fitted distributions per property
└───────────────────────────┘
        │
        ▼
┌───────────────────────────┐
│ 04 generate_monte_carlo_   │──► correlated property realizations
│    realizations            │
└───────────────────────────┘
        │
        ▼
┌───────────────────────────┐
│ 05 build_stratified_gp_    │──► GP-depth-adjusted, correlation-preserved
│    models /                │    realizations
│    integrate_monte_carlo_  │
│    with_gp                 │
└───────────────────────────┘
        │
        ├──► 06 sim_component_comp / simulate_profile_depths_by_mukey
        │    (component composition + horizon depth/thickness variability)
        │
        ├──► 07 calculate_aws_df (available water storage via ROSETTA/Van Genuchten)
        │
        ▼
┌───────────────────────────┐        ┌───────────────────────────┐
│ 09 simulate_ssurgo_mapunit_ │      │ 09 fetch_solus_percentiles │
│    draws -> prior raster   │      │    -> likelihood raster     │
└──────────────┬─────────────┘      └──────────────┬──────────────┘
               │                                    │
               └────────────────┬───────────────────┘
                                ▼
                  ┌───────────────────────────┐
                  │ 09 fuse_property_adaptive  │──► fused posterior SpatRaster
                  │    (built on 08's fusion    │
                  │    primitives, raster-native)│
                  └───────────────────────────┘
                                │
                                ▼
                  ┌───────────────────────────┐
                  │ 02 validate_complete_       │──► quality/validation report
                  │    workflow                 │
                  └───────────────────────────┘
```

Not every path runs through all of this - e.g. a caller who only wants tabular Monte Carlo
estimates from SSURGO alone stops after step 05/06 and never touches the raster fusion group.
`08 Bayesian updating (scalar)` provides the same fusion primitives `09`'s raster pipeline is built
on, but as standalone building blocks, not yet called from `04`'s tabular pipeline directly.

## Group documents

| #  | Document                                                                              | Functional area                                                 | Key R/ files                                                                                                                |
| -- | ------------------------------------------------------------------------------------- | --------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| 01 | [01_data_acquisition_processing.md](01_data_acquisition_processing.md)                 | SSURGO data acquisition, cleaning, infilling                    | `ssurgo-acquisition.R`, `ssurgo-processing.R`, `data-infilling.R`                                                     |
| 02 | [02_statistics_diagnostics.md](02_statistics_diagnostics.md)                           | Statistical characterization + workflow QA                      | `statistics.R`, `validation-diagnostics.R`                                                                              |
| 03 | [03_distribution_fitting_correlations.md](03_distribution_fitting_correlations.md)     | Percentile fitting, ILR transforms, correlation matrices (core) | `distributions.R`, `percentile-sampling.R`, `kssl-reference-correlations.R`                                           |
| 04 | [04_monte_carlo_simulation.md](04_monte_carlo_simulation.md)                           | Correlated Monte Carlo simulation engine (core)                 | `monte-carlo.R`                                                                                                           |
| 05 | [05_gp_modeling_multivariate_adjustment.md](05_gp_modeling_multivariate_adjustment.md) | GP depth-trend modeling + MC/GP integration (core)              | `gp-modeling.R`, `multivariate-adjustment.R`                                                                            |
| 06 | [06_profile_component_depth_simulation.md](06_profile_component_depth_simulation.md)   | Component composition + horizon depth/thickness simulation      | `property-simulation.R`, `depth-simulation.R`                                                                           |
| 07 | [07_aws_van_genuchten_modeling.md](07_aws_van_genuchten_modeling.md)                   | Available water storage / Van Genuchten modeling                | `aws-simulation.R`                                                                                                        |
| 08 | [08_bayesian_updating.md](08_bayesian_updating.md)                                     | Scalar Bayesian updating/fusion (core)                          | `bayesian-updating.R`                                                                                                     |
| 09 | [09_multi_source_raster_fusion_pipeline.md](09_multi_source_raster_fusion_pipeline.md) | Raster-native fusion core + SSURGO/SOLUS100 adapters            | `distribution-fitting-raster.R`, `raster-fusion.R`, `raster-cache.R`, `ssurgo-simulation.R`, `solus-simulation.R` |
| 10 | [10_utilities.md](10_utilities.md)                                                     | Shared validation, logging, config, I/O (leaf module)           | `utils.R`                                                                                                                 |

## Known issues surfaced while writing this documentation

Producing these docs required reading every function body (not just signatures), which turned up a
few real, previously-undocumented issues. Three were fixed (with a regression test and a clean
`devtools::check()` afterward); the fourth turned out not to be a bug on closer inspection:

- **Fixed**: `validate_monte_carlo_quality()` (`validation-diagnostics.R`, see
  `02_statistics_diagnostics.md`) called `generate_simulation_diagnostics()` with 2 arguments, but
  that function (defined in `monte-carlo.R`) requires 5. The call always errored and was silently
  absorbed by a surrounding `tryCatch` into a `validation_failed` placeholder result. Fixed by
  passing all 5 args, sourced from data already in scope (`original_data`, the local
  `simulation_data`, and `monte_carlo_results$metadata$properties`/`$correlation_structure`).
- **Fixed**: `raster-fusion.R`'s top-of-file comment used to say `run_stage1_fusion()`/
  `run_stage1_fusion_group()` were "not ported" because their dependencies "don't exist anywhere in
  this repo" - stale, since all of those dependencies (`raster-cache.R`, `ssurgo-simulation.R`,
  `solus-simulation.R`) now exist and both orchestrator functions are fully implemented later in the
  same file. Updated the comment to describe the current state (see
  `09_multi_source_raster_fusion_pipeline.md`).
- **Fixed**: the `sim_comppct` integration gap between `property-simulation.R` and
  `depth-simulation.R` (see `06_profile_component_depth_simulation.md`) - `simulate_profile_depths_by_mukey()`
  now calls `sim_component_comp()` and left-joins the result onto every horizon row by `cokey`
  internally, mirroring the already-verified pattern in `ssurgo-simulation.R`'s
  `simulate_ssurgo_mapunit_draws()`. Its previously-unused `n_simulations` parameter now flows into
  that call. Verified with a new deterministic offline regression test (mocked SDA fetch, a
  degenerate `comppct_l = comppct_r = comppct_h = 100` component that makes `sim_comppct` exactly
  equal to `n_simulations`). `simulate_and_perturb_soil_profiles()`'s own contract - it still
  requires the join to already be done when called directly - is unchanged.
- **Not a bug, on closer inspection**: `utils.R`'s `get_predefined_properties("ssurgo")` calling
  `create_ssurgo_property_lookup_working()` (which lives in `ssurgo-acquisition.R`) is a real
  exception to "utils.R is a leaf module," but it's a deliberate, already-`tryCatch`-guarded design
  (see the inline comment at that call site) - R packages have no per-file circular-dependency
  restriction within one namespace, so this isn't a functional defect, just an architectural note
  worth knowing about (see `10_utilities.md`). Left as-is, along with `utils.R`'s few unused exported
  functions (`check_required_columns()`, `export_workflow_metadata()`, most of the WKT
  geometry-validation helpers) - unused surface area, not broken code, and removing public API
  without confirming nothing external depends on it isn't a call to make silently.
