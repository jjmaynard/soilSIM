# Build a KSSL-derived fallback correlation matrix for a set of properties

Starts from an identity matrix dimnamed by `properties`, then overlays
the *entire* principal submatrix of the matched KSSL property matrix at
the positions of whichever `properties` have a mapping in
[`.kssl_property_name_map()`](https://jjmaynard.github.io/soilSIM/reference/dot-kssl_property_name_map.md).
Properties without a mapping stay at their identity values (1 on the
diagonal, 0 cross-correlation with everything, including other mapped
properties) - a deliberate, documented partial degrade rather than an
all-or-nothing failure.

## Usage

``` r
build_kssl_fallback_matrix(properties, genhz = NULL)
```

## Arguments

- properties:

  Character vector of property names to build a matrix for (in the order
  the returned matrix should be dimnamed).

- genhz:

  Optional single genhz value (`"O"`/`"A"`/`"E"`/`"B"`/`"C"`/`"Cr"` for
  the property matrix - `"R"` has no property-matrix entry). `NULL`
  (default) pools an unweighted average across every available genhz key
  (positive-definite matrices form a convex cone, so the average stays
  positive-definite).

## Value

A `length(properties) x length(properties)` correlation matrix dimnamed
by `properties`, or `NULL` if fewer than 2 of `properties` have a KSSL
mapping, or if `genhz` is supplied but not among the available keys.

## Details

The result is positive-definite by construction: principal submatrices
of a positive-definite matrix are positive-definite, and a
block-diagonal matrix (KSSL block + identity block, zero cross-terms)
with positive-definite blocks is itself positive-definite. It is still
passed through
[`ensure_positive_definite_matrix()`](https://jjmaynard.github.io/soilSIM/reference/ensure_positive_definite_matrix.md)
defensively (cheap, and
[`configure_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/configure_correlation_structure.md)
already does this again on the final matrix regardless).
