# Calculate Horizon Quantiles and Probabilities by Mukey

Calculates quantiles (5th, 50th, 95th percentiles) and prediction
interval widths for soil properties by mukey and depth, across a wide
property basket (sand/silt/clay, bulk density, water retention, RFV, pH,
CEC, SOC, optional van Genuchten params), and, when texture data and the
optional `soiltexture` package are both available, the most probable
USDA soil texture class. Ported from `code/sim-functions.R:2728`.

## Usage

``` r
hz_quant_prob_mukey(hz_data)
```

## Arguments

- hz_data:

  A data frame of simulated horizon data, with `mukey`, `hzdept_r`,
  `hzdepb_r`, and any of the recognized property columns (see the
  `prop_mapping`/`optional_props` translation table in the source).

## Value

A data frame with one row per `mukey`/depth, quantile statistics
(`_05`/`_50`/`_95`/ `_PIW90` suffixed columns), and - when texture data
and `soiltexture` are available - the most probable texture class and
its simulation-frequency probability.
