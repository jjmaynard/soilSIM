# Soil Profile Depth Simulation

SSURGO-based soil-profile horizon-depth simulation: top-down/bottom-up
depth simulation, thickness-variability estimation via repeated
simulation, OSD-distinctness-derived boundary perturbation, and
[`aqp::SoilProfileCollection`](https://ncss-tech.github.io/aqp/reference/SoilProfileCollection-class.html)-level
orchestration (single profile, whole collection, parallelized
collection, or by-mukey). Ported from
`code_ref/brdf/depth_simulation.R`, the canonical, generalized version
of this logic (superseding earlier copies kept for historical reference
in `code/scratch/`, not part of this port).

This is the first place `aqp`/`SoilProfileCollection` objects are used
in `soilSIM`; the random triangular-distribution sampler these functions
rely on,
[`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md),
lives in `R/distributions.R` (it was originally defined in the
*companion* legacy file, `code_ref/brdf/property_simulation.R`, not in
`depth_simulation.R` itself). That companion file's correlated-property
and van Genuchten/AWS modeling functions remain out of scope for this
port.
