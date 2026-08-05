# Ensure Essential Columns Are Present and Typed

Guarantees `cokey` and `hzname` exist (synthesizing placeholder values
if missing) and coerces `hzdept_r`/`hzdepb_r` to numeric when present.

## Usage

``` r
ensure_essential_columns_working(df)
```

## Arguments

- df:

  Horizon data frame.

## Value

`df` with the essential columns present/coerced.
