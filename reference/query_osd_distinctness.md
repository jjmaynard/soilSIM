# Query OSD Data and Convert Horizon Distinctness to Offset

Retrieves Official Series Description (OSD) data for the soil series
(component names) present in `horizon_data`, extracts horizon
distinctness, and converts the distinctness codes into boundary-offset
values via
[`aqp::hzDistinctnessCodeToOffset()`](https://ncss-tech.github.io/aqp/reference/hzDistinctnessCodeToOffset.html).

## Usage

``` r
query_osd_distinctness(horizon_data)
```

## Arguments

- horizon_data:

  Data frame containing SSURGO horizon data, with at least `compname`
  and `hzname` columns.

## Value

A data frame containing the following columns: - `id`: Unique identifier
for the soil profile (the component name). - `hzname`: The name of the
soil horizon. - `distinctness`: The distinctness code for the horizon
(e.g., A, C, G, D). - `genhz`: The generalized horizon designation. -
`bound_sd`: The calculated offset value from the distinctness code.

## Examples

``` r
if (FALSE) { # \dontrun{
  horizon_data <- data.frame(
    compname = c("amador", "pentz", "pardee", "auburn", "loafercreek", "millvilla"),
    hzname = c("A", "Bt", "R", "A", "Bw", "C")
  )
  result <- query_osd_distinctness(horizon_data)
  head(result)
} # }
```
