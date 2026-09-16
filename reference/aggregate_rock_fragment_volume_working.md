# Aggregate Rock Fragment Volume (Working Version)

Sums `chfrags`' per-fragment-size-class `rfv_l`/`rfv_r`/`rfv_h` values
into one horizon-level (`chkey`) total. A chorizon with zero `chfrags`
rows (genuinely 0% rock fragments - common for fine-textured horizons)
is explicitly set to `rfv_l = rfv_r = rfv_h = 0`, distinguishing it from
a chorizon whose fragment data was never queried at all (which stays
`NA`) - both would otherwise collapse to the same `NA` after the
aggregate's `left_join()`, and downstream
[`impute_rfv_values()`](https://jjmaynard.github.io/soilSIM/reference/impute_rfv_values.md)
treats any `NA` as "missing," fabricating a non-zero estimate for what
was actually a real zero.

## Usage

``` r
aggregate_rock_fragment_volume_working(ssurgo_data, verbose = FALSE)
```

## Arguments

- ssurgo_data:

  SSURGO data frame

- verbose:

  Logical; provide progress messages

## Value

Data frame with aggregated RFV values (0, not `NA`, for genuinely
fragment-free horizons)
