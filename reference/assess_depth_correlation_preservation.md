# Assess Depth-Specific Correlation Preservation

Compares a depth-specific correlation matrix (`depth_cor`) against the
overall/original correlation structure, reusing
[`assess_correlation_differences()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_differences.md)
on the common properties.

## Usage

``` r
assess_depth_correlation_preservation(
  depth_cor,
  original_correlations,
  depth,
  criteria
)
```

## Arguments

- depth_cor:

  Correlation matrix at one depth.

- original_correlations:

  List with a `global_correlation_matrix` entry.

- depth:

  Depth value (used only for logging/context).

- criteria:

  Assessment criteria (passed through to
  [`assess_correlation_differences()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_differences.md)).

## Value

List with `preservation_quality` (1 - mean difference, clamped to
`[0, 1]`), `correlation_differences` (mean absolute difference),
`depth_specific_assessment`.
