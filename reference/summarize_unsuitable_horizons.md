# Summarize Unsuitable Horizons

Restores the reporting half of the legacy `filter_unsuitable_horizons()`
(its flagging half is already handled by
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)):
logs and returns a summary of which horizon types were excluded from
infilling. Uses
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md)
(this package's established logging convention) rather than the legacy
[`cat()`](https://rdrr.io/r/base/cat.html), a deliberate improvement
consistent with how every other diagnostic message in this package is
emitted.

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
