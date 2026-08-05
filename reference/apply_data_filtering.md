# Apply Data Filtering

Applies optional data filtering based on user preferences.

## Usage

``` r
apply_data_filtering(
  df,
  remove_unsuitable,
  remove_incomplete,
  required_properties,
  all_properties,
  verbose
)
```

## Arguments

- df:

  Input data frame

- remove_unsuitable:

  Remove unsuitable horizons

- remove_incomplete:

  Remove incomplete rows

- required_properties:

  Required complete properties

- all_properties:

  All processed properties

- verbose:

  Verbose output

## Value

Filtered data frame
