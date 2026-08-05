# Infill Missing Distinctness Values for Horizons

Fills in missing `distinctness` values for common horizon names (O, A,
B, C, R, Cr) based on typical boundary characteristics observed in the
field.

## Usage

``` r
infill_missing_distinctness(horizon_data)
```

## Arguments

- horizon_data:

  A data frame containing horizon data. It must include `hzname` and
  `distinctness` columns.

## Value

A data frame with missing `distinctness` values infilled based on the
horizon name or generalized horizon group.

## Details

### Default Boundary Distinctness Assignments:

- **O Horizons (Organic Layers):** diffuse

- **A Horizons (Topsoil):** clear

- **B Horizons (Subsoil):** gradual

- **C Horizons (Parent Material):** gradual

- **Cr Horizons (Weathered Bedrock):** gradual

- **R Horizons (Bedrock):** abrupt

## Examples

``` r
if (FALSE) { # \dontrun{
  df <- data.frame(
    hzname = c("A", "B", "C", "R", "O", "Cr"),
    distinctness = c(NA, "gradual", NA, NA, "diffuse", NA),
    stringsAsFactors = FALSE
  )
  df <- infill_missing_distinctness(df)
  print(df)
} # }
```
