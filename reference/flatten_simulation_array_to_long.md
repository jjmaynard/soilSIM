# Flatten a Monte Carlo Simulation Array to Long Format

Converts the `[horizon, property, realization]` array returned by
[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)
(as `result$simulation_data`) into a long-format data frame with one row
per horizon-realization combination, matching the shape
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md),
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md),
and
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md)
expect (a `simulation_number` column plus one column per property,
alongside the original horizon metadata).

## Usage

``` r
flatten_simulation_array_to_long(sim_array, cokey_data)
```

## Arguments

- sim_array:

  A `[horizon, property, realization]` array as produced by
  [`simulate_correlated_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_properties.md).

- cokey_data:

  The original per-horizon data frame passed in as `soil_data` - must
  have exactly one row per `sim_array` horizon index, in the same order.

## Value

Long-format data frame.
