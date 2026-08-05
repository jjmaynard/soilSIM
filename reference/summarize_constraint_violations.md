# Summarize Constraint Violations

Aggregates the nested property -\> group -\>
[`assess_trend_realism()`](https://jjmaynard.github.io/soilSIM/reference/assess_trend_realism.md)
results into an overall violation count and compliance rate.

## Usage

``` r
summarize_constraint_violations(trend_predictions)
```

## Arguments

- trend_predictions:

  Nested list: property -\> group -\> trend assessment (as returned by
  [`assess_trend_realism()`](https://jjmaynard.github.io/soilSIM/reference/assess_trend_realism.md)).

## Value

List with `total_violations`, `violation_types` (properties with at
least one non-realistic trend), `overall_compliance`.
