# Remove Invalid Component Records

Drops component rows with a missing/empty `cokey` or an out-of-range
`comppct_r` (component percentage, expected `[0, 100]`).

## Usage

``` r
remove_invalid_components_working(df, verbose = FALSE)
```

## Arguments

- df:

  Component data frame.

- verbose:

  Logical; log the number of rows removed.

## Value

`df` with invalid rows removed.
