# KSSL texture (sand/silt/clay) correlation matrices (internal accessor)

Unused by
[`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)
in this version - `ilr1`/`ilr2` are sourced from
[`.kssl_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/dot-kssl_property_matrices.md)'s
bundled columns instead, since that matrix has solid eigenvalues
(0.08-0.22) across every genhz group, unlike this standalone 3x3 matrix
(exactly singular for `E`/`Cr`/`R`, marginally positive-definite for the
rest). Kept available for any future use.

## Usage

``` r
.kssl_texture_matrices()
```

## Value

Named list, one 3x3 correlation matrix per genhz key
(`O`,`A`,`E`,`B`,`C`,`Cr`,`R`).
