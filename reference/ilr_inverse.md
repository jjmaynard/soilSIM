# Inverse ILR transform: 2 unconstrained coordinates -\> a valid composition

Output columns `clay`/`sand`/`silt` are positions 1/2/3 of the same
sequential binary partition
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)
uses - see that function's doc and the file header comment above for the
positional-role vs identity distinction.

## Usage

``` r
ilr_inverse(z1, z2, total = 100)
```

## Arguments

- z1, z2:

  Numeric vectors of ILR coordinates.

- total:

  The composition's target sum (default 100).

## Value

A 3-column matrix (`clay`, `sand`, `silt`) summing to `total`.
