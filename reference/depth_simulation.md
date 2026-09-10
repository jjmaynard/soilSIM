# Soil Profile Depth Simulation

SSURGO-based soil-profile horizon-depth simulation: top-down/bottom-up
depth simulation, thickness-variability estimation via repeated
simulation, OSD-distinctness-derived boundary perturbation, and
[`aqp::SoilProfileCollection`](https://ncss-tech.github.io/aqp/reference/SoilProfileCollection-class.html)-level
orchestration (single profile, whole collection, parallelized
collection, or by-mukey).

This is the first place `aqp`/`SoilProfileCollection` objects are used
in `soilSIM`. The random triangular-distribution sampler these functions
rely on,
[`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md),
lives in `R/distributions.R`.
