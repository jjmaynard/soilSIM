# Remove Organic Layers and Adjust Depths

Removes rows with organic horizons (where `hzname` contains a capital
"O") and recalculates the remaining horizons' depths within each `cokey`
group so they stay cumulative and re-anchored to 0, preserving the
original horizon thicknesses.

## Usage

``` r
remove_organic_layer(df)
```

## Arguments

- df:

  A data frame containing soil horizon data, with columns `cokey`,
  `hzname`, `hzdept_r`, `hzdepb_r`, and optionally `hzdept_l`,
  `hzdepb_l`, `hzdept_h`, `hzdepb_h`.

## Value

A data frame with organic layers removed and depths adjusted within each
`cokey` group.
