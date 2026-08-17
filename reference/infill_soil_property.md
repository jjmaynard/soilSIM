# Main Soil Property Infilling Function

Enhanced version of the core infilling function with comprehensive
recovery strategies and automatic exclusion of unsuitable horizons (R,
Cr, O horizons).

## Usage

``` r
infill_soil_property(
  df,
  property_name,
  property_config = NULL,
  max_depth = 250,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- df:

  Input soil data frame

- property_name:

  Name of the property to infill

- property_config:

  Optional property configuration

- max_depth:

  Maximum depth for infilling (default: 250 cm)

- verbose:

  Whether to provide detailed progress messages

## Value

Data frame with infilled property data
