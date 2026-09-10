# Map soilSIM property names to KSSL reference-matrix column names

`om` -\> `soc`: the KSSL matrix's `soc` column was fit on organic-matter
(`om_r`) data, not true lab-measured soil organic carbon - so it's an
organic-matter *correlation structure*, not a true-SOC one. The
simulated `soc` values are not on the OM scale:
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
rescales the `om` triplet to an estimated SOC triplet
(`OM_TO_SOC_FACTOR`, the Van Bemmelen factor) before it enters the
correlated draw - the matrix itself needed no change, since Pearson
correlation is invariant under a positive linear rescale of one
variable. Treat the *correlations* involving `soc` as an
organic-matter-based approximation (real SOC correlations may differ
slightly), but treat the simulated `soc` *values* as genuine SOC-scale
estimates, not raw OM.

## Usage

``` r
.kssl_property_name_map
```

## Format

An object of class `character` of length 14.

## Details

`ilr1`/`ilr2` map directly (same names) - valid only because
`monte-carlo.R`'s `composition_groups$texture$members` default is
`(sandtotal, silttotal, claytotal)`, matching the sequential binary
partition `compositions::ilr()` uses to build the KSSL matrix (sand vs
silt+clay, then silt vs clay). See `distributions.R`'s ILR section
header for the positional-role convention this depends on.

`caco3`/`ec`/`ecec`/`gypsum`/`sar` -\> themselves: these 5 chemistry
properties have **no** entry in
[`.kssl_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/dot-kssl_property_matrices.md) -
the identity passthrough exists so
[`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)
recognizes them as "mapped" property names, but its own
`mapped_kssl_names %in% rownames(kssl_source)` check (see that
function's implementation) then finds no matching row in the actual KSSL
matrix and silently excludes them - they get
[`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)'s
identity-matrix default (uncorrelated with everything, including each
other) rather than an error. **Upgrade path**: once raw KSSL lab data
covering these 5 properties is available, re-run
`data-raw/build_kssl_reference_correlations.R` to produce an expanded
`kssl_property_matrices` sysdata object whose rows/columns include
them -
[`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)
picks the real correlations up automatically via that same check, with
**no code changes** needed here or anywhere else that calls it. Keep any
future real-matrix column names identical to these 5 (plain property
name, no rename) so the mapping below doesn't need to change too.
