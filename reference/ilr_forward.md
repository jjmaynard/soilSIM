# Forward ILR transform: 2 unconstrained coordinates from a 3-part composition

`clay`/`sand`/`silt` name positions 1/2/3 of the sequential binary
partition (position 1 vs positions 2 and 3, then 2 vs 3) - which real
property occupies which position is entirely up to the caller (see the
file header comment above); the names are historical, not an identity
requirement.

## Usage

``` r
ilr_forward(clay, sand, silt)
```

## Arguments

- clay, sand, silt:

  Numeric vectors, same units (e.g. percent), same length.

## Value

A 2-column matrix (`z1`, `z2`).
