# Remove Invalid Horizon Records

Drops horizon rows with missing/inconsistent depths (negative top depth,
bottom not greater than top, either depth beyond `max_depth`) or a
missing/empty `cokey`.

## Usage

``` r
remove_invalid_horizons_working_compatible(
  df,
  max_depth = 250,
  verbose = FALSE
)
```

## Arguments

- df:

  Horizon data frame.

- max_depth:

  Maximum plausible depth (cm).

- verbose:

  Logical; log the number of rows removed.

## Value

`df` with invalid rows removed.
