# Water Retention Estimation using Saxton-Rawls

Estimates missing water retention values using Saxton-Rawls equations
with rock fragment correction, integrated with the main workflow.

## Usage

``` r
infill_water_retention_saxton_rawls_integrated(
  df,
  max_depth = 250,
  add_ranges = TRUE,
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- df:

  Input data frame

- max_depth:

  Maximum depth for estimation

- add_ranges:

  Whether to add range estimates

- overwrite:

  Whether to overwrite existing values

- verbose:

  Whether to provide progress messages

## Value

Data frame with estimated water retention values
