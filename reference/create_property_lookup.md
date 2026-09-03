# Simple Property Lookup Creator

Helper function to create custom property lookups easily

## Usage

``` r
create_property_lookup(properties, synonyms = NULL, metadata = NULL)
```

## Arguments

- properties:

  Vector of property names

- synonyms:

  Optional named list of synonyms

- metadata:

  Optional metadata list

## Value

Property lookup object

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API
for callers who need it, not dead/broken code.
