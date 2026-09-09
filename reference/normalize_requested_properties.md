# Normalize a Requested-Properties Vector to `simulate_cokey_generalized()` Parameter Names

Resolves a caller's mixed vocabulary (caller-facing ids like
`"clay"`/`"bulk_density"`, SSURGO stems like `"ph1to1h2o"`,
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
output columns like `"clay_total"`, or the internal parameter names
themselves) into the fixed
`c("db","wr_3b","wr_15b","ilr1","ilr2","rfv","ph","cec","soc")`
vocabulary that
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
gates its per-row `param_list` on.

## Usage

``` r
normalize_requested_properties(requested_properties)
```

## Arguments

- requested_properties:

  `NULL`, or a character vector of property identifiers in any of the
  accepted vocabularies.

## Value

`NULL`, or a de-duplicated character vector in canonical `param_order`
order.

## Details

Texture is a coupled all-or-nothing unit: any texture member
(`sand`/`silt`/`clay` in any spelling, or `ilr1`/`ilr2` directly) pulls
in **both** ILR axes, since
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
derives `ilr1`/`ilr2` jointly from one sand+silt+clay draw and only ever
adds them as a pair. The only
intra-[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
derivation is this texture coupling - any future dependency (e.g. a
property estimated from another) is expanded here.

Total and idempotent:
`normalize_requested_properties(normalize_requested_properties(x))`
equals `normalize_requested_properties(x)`, so callers (including
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
itself) can normalize freely. `NULL` in -\> `NULL` out ("simulate
everything", the unchanged default).
