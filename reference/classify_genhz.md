# Classify a horizon name into a generalized master-horizon designation

Strips a leading lithologic-discontinuity digit prefix (e.g. `"2Bt2"`
-\> `"Bt2"`), then matches against the master-horizon letters the KSSL
reference correlation matrices are keyed by: `O`, `A`, `E`, `B`, `C`,
`Cr`, `R`. `Cr` is checked before `C` so `"Cr"`/`"Crt"` aren't
misclassified as `"C"`.

## Usage

``` r
classify_genhz(hzname)
```

## Arguments

- hzname:

  Character vector of raw SSURGO horizon-designation text (e.g. `"Bt2"`,
  `"2Bw"`, `"Cr"`, `"R"`). `NA` is preserved as `NA`.

## Value

Character vector, same length as `hzname`, of
`O`/`A`/`E`/`B`/`C`/`Cr`/`R` or `NA_character_` for unrecognized text.
