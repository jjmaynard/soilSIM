# Assess Simulation Quality

Derives real constraint-satisfaction and statistical-consistency scores
from `diagnostics$per_property` (simulated vs. original mean/SD), and
combines them with `output_validation$success_rate` into an overall
quality score.

## Usage

``` r
assess_simulation_quality(
  simulation_results,
  diagnostics,
  output_validation,
  config
)
```

## Arguments

- simulation_results:

  Unused (kept for interface compatibility).

- diagnostics:

  Result of
  [`generate_simulation_diagnostics()`](https://jjmaynard.github.io/soilSIM/reference/generate_simulation_diagnostics.md).

- output_validation:

  Result of
  [`validate_simulation_output()`](https://jjmaynard.github.io/soilSIM/reference/validate_simulation_output.md),
  with a `success_rate` field.

- config:

  Unused (kept for interface compatibility).

## Value

List with `overall_quality_score`, `component_scores`
(`data_quality`/`constraint_satisfaction`/`statistical_consistency`),
`recommendations`.
