# Saxton-Rawls water retention (field capacity, wilting point), raster / stack-native

The vectorized, `terra`-native counterpart of
[`calculate_saxton_rawls_single()`](https://jjmaynard.github.io/soilSIM/reference/calculate_saxton_rawls_single.md) -
identical equations and clamps, but every operation is `terra`
`Arith`/`Math`/`clamp`/`ifel` so it runs on multi-layer realization
stacks in one pass. Used by
[`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)'s
default (`method = "saxton_rawls"`) path, because SOLUS100 publishes
**no** water-retention variable - `wr_3b`/`wr_15b` can never be fused
directly, only derived from the fusable
sand/silt/clay/`dbovendry`/`soc`/`fragvol`.

## Usage

``` r
saxton_rawls_raster(sand, clay, silt, db, rfv, om)
```

## Arguments

- sand, clay, silt, db, rfv, om:

  Percent (or g/cm^3 for `db`) `SpatRaster`s. Any may be single or
  multi-layer; single-layer inputs are recycled.

## Value

`list(fc =, wp =)` - volumetric % `SpatRaster`s, `wp < fc` enforced.

## Details

Because
[`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md)
re-marginalizes `sand`/`silt`/`clay` to their own fused posteriors
independently, a realization's texture triple may not sum to 100. This
reproduces
[`calculate_saxton_rawls_single()`](https://jjmaynard.github.io/soilSIM/reference/calculate_saxton_rawls_single.md)'s
renormalization: where the texture sum is off by more than 5 points, the
three fractions are rescaled to sum to 100 before the equations run (the
equations themselves use only the sand and clay fractions).
