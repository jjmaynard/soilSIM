# Infer Simulation Property Columns

Infers which columns of a long-format simulation data frame represent
soil properties to adjust, as the numeric columns remaining after
excluding known structural/metadata columns.

## Usage

``` r
infer_simulation_properties(simulation_data)
```

## Arguments

- simulation_data:

  Long-format simulation data frame.

## Value

Character vector of inferred property column names.
