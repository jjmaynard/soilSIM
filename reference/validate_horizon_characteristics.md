# Validate Horizon Characteristics

Checks whether surface (`hzdept_r <= 30`) and subsurface horizons fall
within realistic property ranges (`criteria$realistic_ranges`,
defaulting to `get_realistic_property_ranges()`), and measures the
smoothness of the depth transition for the first available property.

## Usage

``` r
validate_horizon_characteristics(simulation_data, criteria)
```

## Arguments

- simulation_data:

  Simulation data.

- criteria:

  List, optionally with `realistic_ranges`.

## Value

List with `surface_horizon_quality`, `subsurface_quality`,
`transition_quality` (each a fraction/score in `[0, 1]`, or `NA` when
not computable).
