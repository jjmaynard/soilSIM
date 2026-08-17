# Calculate Derived Horizon Properties

Fills in `hzthk_r` (thickness), `hz_midpoint`, and `awc_r` (available
water capacity, `wthirdbar_r - wfifteenbar_r`) when the source columns
they're derived from are present and the derived column is not already
present.

## Usage

``` r
calculate_derived_horizon_properties_working(df, verbose = FALSE)
```

## Arguments

- df:

  Horizon data frame.

- verbose:

  Logical; log a completion message.

## Value

`df` with the derived columns added.
