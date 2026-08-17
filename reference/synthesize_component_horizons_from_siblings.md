# Synthesize a Representative Horizon Profile for a Component from AOI Siblings

Given one component known to exist (via a real `comppct`) but missing
all horizon data (`target_component`), and the already fully processed
horizon rows of OTHER cokeys within the same
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
call's `ssurgo_data` that share `target_component$compname`
(`sibling_horizons`), builds a synthetic horizon profile for
`target_component$cokey` by averaging sibling horizons within each
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
group.

## Usage

``` r
synthesize_component_horizons_from_siblings(target_component, sibling_horizons)
```

## Arguments

- target_component:

  One-row data frame: `mukey`, `cokey`, `compname`, `comppct_l`,
  `comppct_r`, `comppct_h`.

- sibling_horizons:

  Data frame of horizon rows (same column shape as
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
  `ssurgo_data`) from one or more OTHER cokeys sharing
  `target_component$compname`.

## Value

Data frame of synthesized horizon rows for `target_component$cokey`, one
row per genhz group produced, with
`infill_method = "component_aoi_average"` and
`component_synthesized = TRUE` on every row. A 0-row data frame if no
sibling horizon was classifiable at all.

## Alignment method

Horizons are aligned by `classify_genhz(hzname)` group
(`O`/`A`/`E`/`B`/`C`/`Cr`/`R`), not by position or depth bin - robust to
sibling profiles having different horizon counts or depths. Sibling rows
whose `hzname` doesn't classify
([`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
returns `NA`) are dropped from averaging entirely rather than guessed
at - verified before choosing this: no depth-or-property- based genhz
fallback exists anywhere else in soilSIM to reuse
([`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
itself and the separate
[`aqp::generalizeHz()`](https://ncss-tech.github.io/aqp/reference/generalize.hz.html)
usage in `R/depth-simulation.R` are both purely `hzname`-regex
matchers).

## Averaging rule per genhz group

- Numeric columns (depths, thicknesses, every property column present):
  arithmetic mean across whichever sibling rows fall in that genhz
  group, `na.rm = TRUE`.

- `hzname`: the most frequent `hzname` string among the group's sibling
  rows (ties broken alphabetically) - reused verbatim, never invented,
  so it re-classifies to the same genhz group deterministically.

- Other character/logical columns (`texcl`, `desgnmaster`,
  `unsuitable_horizon`, etc.): most frequent non-`NA` value in the
  group; `NA` if every sibling value in the group is `NA`.

- `mukey`/`cokey`/`compname`/`comppct_l`/`comppct_r`/`comppct_h`: taken
  from `target_component`, never averaged from siblings - the missing
  component's own real identity and composition percentages are already
  known and correct.

- `chkey`: freshly generated (`"synth_<cokey>_<genhz>"`), unique and
  clearly synthetic.

A genhz group present in only one sibling still produces a row (mean of
one = a copy); a group present in zero siblings is simply absent from
the output - never fabricated from nothing.
