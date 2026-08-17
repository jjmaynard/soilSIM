# Execute SSURGO Query (Working Version)

Executes SSURGO database query exactly as in working code

## Usage

``` r
execute_ssurgo_query_working(
  mukey_list,
  properties,
  ssurgo_lookup,
  include_restrictions = TRUE,
  verbose = FALSE
)
```

## Arguments

- mukey_list:

  Vector of map unit keys

- properties:

  Vector of properties to retrieve

- ssurgo_lookup:

  Property lookup table

- include_restrictions:

  Logical; include restriction data

- verbose:

  Logical; provide progress messages

## Value

Data frame with SSURGO data
