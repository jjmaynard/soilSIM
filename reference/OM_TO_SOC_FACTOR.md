# Van Bemmelen-factor conversion from SSURGO organic matter % to an estimated soil organic carbon %

SSURGO's `chorizon` table reports organic matter (`om_l/_r/_h`), not
organic carbon - but
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
simulated `soc` column (and the KSSL reference correlation matrix's
`soc` row/column, see `R/kssl-reference-correlations.R`'s
`.kssl_property_name_map` docs) both need an actual SOC-scale value to
be fused against SOLUS100's real `soc` variable (organic carbon).
`SOC = OM / 1.724` is the standard Van Bemmelen-factor approximation
(the same factor
[`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)'s
`soc_to_om` argument uses in the opposite direction,
`R/raster-fusion-bridge.R`). Applying this conversion needs no change to
the KSSL correlation *matrix* itself - Pearson correlation is invariant
under a positive linear rescale of one variable, and
triangular-distribution parameters (min/mode/max) transform linearly
too - only the marginal values entering the correlated draw need
rescaling, which is exactly where this constant is used.

## Usage

``` r
OM_TO_SOC_FACTOR
```

## Format

An object of class `numeric` of length 1.
