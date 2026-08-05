# Map a Property Id to its Simulated-Data Column Name

Map a Property Id to its Simulated-Data Column Name

## Usage

``` r
property_to_sim_column(property_id)
```

## Arguments

- property_id:

  One of `"ph"`, `"bulk_density"`, `"soc"`, `"cec"`, `"clay"`, `"sand"`,
  `"silt"`, `"rock_fragments"`.

## Value

The corresponding column name in
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
output.
