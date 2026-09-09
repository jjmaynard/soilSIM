# SOLUS100's Right-Censoring Depth for `anylithicdpt`/`resdept` (S2)

Both depth-to-bedrock (`anylithicdpt`) and depth-to-restriction
(`resdept`) are right-censored at 201 cm in SOLUS100's own training
data, per gSSURGO data-definition convention (Nauman et al. 2024,
SOLUS100 documentation section 2.1.3, confirmed against the actual
training-data-prep source). A predicted value at/above this censor point
means "no restriction detected within the training depth range" - not a
literal claim that bedrock/restriction sits at exactly 201 cm.
[`fetch_solus_restriction_depth()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_restriction_depth.md)
recodes values `>= SOLUS_RESTRICTION_CENSOR_CM` to `Inf` (never
truncates) rather than treating 201 as a real cutoff.

## Usage

``` r
SOLUS_RESTRICTION_CENSOR_CM
```

## Format

An object of class `numeric` of length 1.
