# Van Genuchten / ROSETTA-Based Available Water Storage (AWS) Modeling

Closed-form van Genuchten water-retention evaluation, Monte Carlo
simulation of available water holding capacity (AWHC) from
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)-derived
pedotransfer parameters, and depth-sliced available-water-storage
summaries for a set of soil components. Ported from
`code_ref/brdf/property_simulation.R` (identical to the superseded
`code/sim-functions.R` copy). Self-contained: unlike
`R/depth-simulation.R`'s functions, none of these three depend on any
other not-yet-ported helper from that file (e.g.
[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)).
