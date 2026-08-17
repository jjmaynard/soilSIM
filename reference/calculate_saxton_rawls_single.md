# Saxton-Rawls Water Retention Pedotransfer Function (Single Horizon)

Estimates field capacity and wilting point water content from texture,
bulk density, rock fragment volume, and organic matter, via the
Saxton-Rawls pedotransfer equations (simplified/RFV-corrected variant).
Real, checkable math (not a placeholder) - inputs are clamped to
physically plausible ranges and texture percentages are renormalized to
sum to 100 when off by more than 5 points.

## Usage

``` r
calculate_saxton_rawls_single(
  sand_pct,
  clay_pct,
  silt_pct,
  bulk_density,
  rfv_pct = 0,
  om_pct = 2
)
```

## Arguments

- sand_pct, clay_pct, silt_pct:

  Texture percentages (0-100).

- bulk_density:

  Bulk density (g/cm^3), clamped to `[0.6, 2.5]`.

- rfv_pct:

  Rock fragment volume percentage, clamped to `[0, 95]`.

- om_pct:

  Organic matter percentage, clamped to `[0.1, 50]`.

## Value

A list of estimated water retention values (field capacity, wilting
point, and their `_l/_h` spread) - see the function body for the exact
returned fields.
