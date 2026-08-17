# Infill Missing Soil Property Values

Loops
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)
(`R/data-infilling.R`, which already special-cases `"rfv"` internally)
over the standard SSURGO property set.

## Usage

``` r
infill_soil_data(df)
```

## Arguments

- df:

  A horizon data frame, as returned by
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md).

## Value

`df` with missing values infilled where possible.
