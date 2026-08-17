# Infill Missing Depth Variability

Infills missing depth variability values in a horizon data frame by
replacing missing top and bottom depth bounds with the representative
depth value plus or minus 2, ensuring that the resulting values are not
negative.

## Usage

``` r
infill_missing_depth_variability(horizon_data)
```

## Arguments

- horizon_data:

  A data frame containing horizon depth data. It must include the
  following columns:

  - `hzdept_r`: Representative top depth.

  - `hzdept_l`: Top low value (may be missing).

  - `hzdept_h`: Top high value (may be missing).

  - `hzdepb_r`: Representative bottom depth.

  - `hzdepb_l`: Bottom low value (may be missing).

  - `hzdepb_h`: Bottom high value (may be missing).

## Value

A data frame with missing depth variability values filled in. Missing
top low values are replaced with `hzdept_r - 2` and missing top high
values with `hzdept_r + 2`. Similarly, missing bottom low values are
replaced with `hzdepb_r - 2` and missing bottom high values with
`hzdepb_r + 2`. All values are ensured to be non-negative.

## Examples

``` r
data <- data.frame(
  hzdept_r = c(10, 20),
  hzdept_l = c(NA, 18),
  hzdept_h = c(NA, NA),
  hzdepb_r = c(30, 40),
  hzdepb_l = c(NA, NA),
  hzdepb_h = c(NA, 42)
)
infill_missing_depth_variability(data)
#>   hzdept_r hzdept_l hzdept_h hzdepb_r hzdepb_l hzdepb_h
#> 1       10        8       12       30       28       32
#> 2       20       18       22       40       38       42
```
