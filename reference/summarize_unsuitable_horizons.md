# Summarize Unsuitable Horizons

Logs and returns a summary of which horizon types were excluded from
infilling (the flagging itself is handled by
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)).

## Usage

``` r
summarize_unsuitable_horizons(df, hzname_col = "hzname")
```

## Arguments

- df:

  Data frame with an `unsuitable_horizon` logical column (see
  [`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md))

- hzname_col:

  Name of the horizon-name column

## Value

List with `n_unsuitable` (count) and `horizon_types` (unique excluded
horizon names)
