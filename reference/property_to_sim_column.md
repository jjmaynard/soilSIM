# Map a Property Id to its Simulated-Data Column Name

Map a Property Id to its Simulated-Data Column Name

## Usage

``` r
property_to_sim_column(property_id)
```

## Arguments

- property_id:

  One of `"ph"`, `"ph1to1h2o"`, `"bulk_density"`, `"dbovendry"`,
  `"soc"`, `"om"`, `"cec"`, `"cec7"`, `"clay"`, `"claytotal"`, `"sand"`,
  `"sandtotal"`, `"silt"`, `"silttotal"`, `"rock_fragments"`, or
  `"rfv"`.

## Value

The corresponding column name in
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
output.
