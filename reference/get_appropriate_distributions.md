# Determine candidate distributions to fit for a property

Computes the existing type-appropriate heuristic candidate list, then,
if `config$distribution_methods` is supplied, narrows to the
intersection of the two - keeping the heuristic as a soft prior rather
than silently discarding an explicit user request. If the intersection
is empty (the user asked only for families the heuristic wouldn't have
suggested), falls back to the user's list unfiltered, with a logged
warning, rather than ignoring it entirely. Previously `config` was
accepted by every caller in this chain but never actually consulted
here - `distribution_methods` was inert.

## Usage

``` r
get_appropriate_distributions(property_name, values, config = NULL)
```

## Arguments

- property_name:

  Property name (drives the type-appropriate heuristic).

- values:

  Numeric vector of observed property values (currently unused by the
  heuristic itself, but taken for interface symmetry with callers and
  potential future data-driven refinement).

- config:

  Optional analysis configuration; reads `config$distribution_methods`.

## Value

Character vector of candidate distribution names.
