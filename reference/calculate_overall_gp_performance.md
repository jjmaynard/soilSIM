# Calculate Overall GP Performance

Aggregates the per-property, per-group results of
[`assess_single_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_gp_performance.md)
into overall mean RMSE/R-squared.

## Usage

``` r
calculate_overall_gp_performance(individual_performance)
```

## Arguments

- individual_performance:

  Nested list: property -\> group -\> performance result (as returned by
  [`assess_single_gp_performance()`](https://jjmaynard.github.io/soilSIM/reference/assess_single_gp_performance.md)).

## Value

List with `mean_r_squared`, `mean_rmse`, `overall_quality`,
`performance_summary`.
