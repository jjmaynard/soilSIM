# Fetch a SOLUS100 Restriction-Depth Raster, Right-Censoring-Guarded (S2)

The truncation-depth signal for bedrock/restriction-aware AWC
([`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)'s
`restriction_depth` argument, MULTI_PROPERTY_FUSION_PLAN.md task S2).
Default variable is `"anylithicdpt"` - trained (per SOLUS100's own
published methodology) **only** on
`reskind %in% c("Lithic bedrock", "Paralithic bedrock")` records, making
it hard-bedrock-specific by construction; no separate SSURGO-side
hard/soft classifier is needed since SOLUS already encodes that
distinction across its two site-level variables. `"resdept"` is trained
on *every* restriction record regardless of kind (fragipans, duripans,
petrocalcic, etc. included) - a broader, more conservative (more
truncation-prone) alternative, offered but not the default, since it
would truncate through partially-penetrable layers with no way to tell
those apart from true bedrock.

## Usage

``` r
fetch_solus_restriction_depth(aoi_vect, variable = "anylithicdpt")
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- variable:

  `"anylithicdpt"` (default, hard-bedrock-specific) or `"resdept"`
  (broader - any restriction, hard or soft).

## Value

A single-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of restriction depth in cm, with every cell
`>= SOLUS_RESTRICTION_CENSOR_CM` recoded to `Inf` (no restriction
detected - never truncates), or `NULL` if the underlying `fetchSOLUS()`
prediction fetch fails.

## See also

`SOLUS_RESTRICTION_CENSOR_CM`,
[`fetch_solus_site_level()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_site_level.md),
[remarginalized_awc()](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)'s
`restriction_depth` argument
