# Apply Composition Constraints

Clamps `normalized_realizations` (an `[n_components, n_realizations]`
matrix) to `constraints$min`/`constraints$max` when supplied.

## Usage

``` r
apply_composition_constraints(normalized_realizations, constraints)
```

## Arguments

- normalized_realizations:

  Matrix of normalized component realizations.

- constraints:

  List, optionally with `min`/`max`, or `NULL`.

## Value

Constrained matrix (same dimensions as input).
