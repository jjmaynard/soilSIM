# Map soilSIM property names to KSSL reference-matrix column names

Confirmed against the exact lookup table that built the source data
(`code_ref/brdf/property_simulation.R:804-816`).

## Usage

``` r
.kssl_property_name_map
```

## Format

An object of class `character` of length 9.

## Details

`om` -\> `soc`: **caveat** - the KSSL matrix's `soc` column is actually
populated from `om_r` (organic matter %) in the source codebase, not
true lab-measured soil organic carbon. Treat it as an organic-matter
proxy, not true SOC.

`ilr1`/`ilr2` map directly (same names) - valid only because
`monte-carlo.R`'s `composition_groups$texture$members` default is
`(sandtotal, silttotal, claytotal)`, matching the sequential binary
partition `compositions::ilr()` used to build the KSSL matrix (sand vs
silt+clay, then silt vs clay). See `distributions.R`'s ILR section
header for the positional-role convention this depends on.
