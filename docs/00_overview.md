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
architecture is split this way so that additional data sources (e.g. HWSD, SoilGrids) can
plug in as new adapters alongside SSURGO/SOLUS, reusing the same simulation/fusion core rather than
duplicating it. See `?soilSIM` for the package-level roxygen summary of this framing.

These 10 documents go one level deeper than that package doc or the per-function
`.Rd`/roxygen reference: each covers a functional group's full function signatures, internal call
graphs, cross-file dependencies, data flow, and known limitations.

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
                      │     fitting-raster.R, core-fusion.R,    │
                      │     cache.R)                       │
                      └───────────────┬───────────────────────────┘
                                      │ consumed by / feeds
              ┌───────────────────────┼───────────────────────┐
              │                       │                       │
   ┌──────────▼──────────┐ ┌──────────▼──────────┐ ┌──────────▼──────────┐
   │  SSURGO ADAPTER      │ │  SOLUS100 ADAPTER    │ │  (further adapters:  │
   │  01 Acquisition/     │ │  09 solus-           │ │   HWSD, SoilGrids,   │
   │     processing        │ │     simulation.R      │ │   etc.)              │
   │  06 Property/depth/   │ │  supplies raster      │ │                      │
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
        ├──► 06 simulate_component_composition / simulate_profile_depths_by_mukey
        │    (component composition + horizon depth/thickness variability)
        │
        ├──► 07 compute_aws (available water storage via ROSETTA/Van Genuchten)
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
on, as standalone building blocks; `04`'s tabular pipeline does not call them directly.

## Group documents

| #  | Document                                                                              | Functional area                                                 | Key R/ files                                                                                                                |
| -- | ------------------------------------------------------------------------------------- | --------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| 01 | [01_data_acquisition_processing.md](01_data_acquisition_processing.md)                 | SSURGO data acquisition, cleaning, infilling                    | `adapter-ssurgo-acquire.R`, `adapter-ssurgo-process.R`, `adapter-ssurgo-infill.R`                                                     |
| 02 | [02_statistics_diagnostics.md](02_statistics_diagnostics.md)                           | Statistical characterization + workflow QA                      | `statistics.R`, `diagnostics.R`                                                                              |
| 03 | [03_distribution_fitting_correlations.md](03_distribution_fitting_correlations.md)     | Percentile fitting, ILR transforms, correlation matrices (core) | `core-distributions.R`, `core-distributions.R`, `core-correlations.R`                                           |
| 04 | [04_monte_carlo_simulation.md](04_monte_carlo_simulation.md)                           | Correlated Monte Carlo simulation engine (core)                 | `core-montecarlo.R`                                                                                                           |
| 05 | [05_gp_modeling_multivariate_adjustment.md](05_gp_modeling_multivariate_adjustment.md) | GP depth-trend modeling + MC/GP integration (core)              | `core-gp.R`, `core-gp.R`                                                                            |
| 06 | [06_profile_component_depth_simulation.md](06_profile_component_depth_simulation.md)   | Component composition + horizon depth/thickness simulation      | `core-simulation.R`, `core-simulation.R`                                                                           |
| 07 | [07_aws_van_genuchten_modeling.md](07_aws_van_genuchten_modeling.md)                   | Available water storage / Van Genuchten modeling                | `model-aws.R`                                                                                                        |
| 08 | [08_bayesian_updating.md](08_bayesian_updating.md)                                     | Scalar Bayesian updating/fusion (core)                          | `core-fusion.R`                                                                                                     |
| 09 | [09_multi_source_raster_fusion_pipeline.md](09_multi_source_raster_fusion_pipeline.md) | Raster-native fusion core + SSURGO/SOLUS100 adapters            | `core-distributions-raster.R`, `core-fusion.R`, `cache.R`, `adapter-ssurgo-simulate.R`, `adapter-solus.R` |
| 10 | [10_utilities.md](10_utilities.md)                                                     | Shared validation, logging, config, I/O (leaf module)           | `utils.R`                                                                                                                 |

