# Validate a `run_stage1_fusion_group()`-shaped cached group result

The
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
counterpart to
[`is_valid_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/is_valid_percentile_list.md) -
see its docs for why an unvalidated cache hit is dangerous. A valid
group result is a named list keyed by every `member_ids` entry, each
element carrying a non-`NULL` `posterior`/`dist` (see
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)'s
return-value docs for the full per-member shape).

## Usage

``` r
is_valid_group_result(x, member_ids)
```

## Arguments

- x:

  A candidate value read from
  [`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
  already
  [`unwrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)-ed.

- member_ids:

  Character vector of expected member ids (from
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)).

## Value

`TRUE` if `x` has the expected shape, else `FALSE`.
