# Create Download Metadata

Creates comprehensive metadata for download operations.

## Usage

``` r
create_download_metadata(
  start_time,
  end_time,
  aoi_wkt,
  properties,
  include_restrictions,
  mukey_count,
  data_rows,
  unique_cokeys,
  spatial_result,
  validation_results,
  components_missing_horizons = NULL,
  components_recovered = NULL
)
```

## Arguments

- start_time:

  Download start time

- end_time:

  Download end time

- aoi_wkt:

  Area of interest WKT

- properties:

  Requested properties

- include_restrictions:

  Whether restrictions were included

- mukey_count:

  Number of map unit keys

- data_rows:

  Number of data rows returned

- unique_cokeys:

  Number of unique component keys

- spatial_result:

  Spatial processing results

- validation_results:

  Data validation results

- components_missing_horizons:

  Data frame or `NULL` - components with a real `comppct` but no AOI
  sibling to recover a profile from (see
  [`recover_missing_horizon_components()`](https://jjmaynard.github.io/soilSIM/reference/recover_missing_horizon_components.md)).

- components_recovered:

  Data frame or `NULL` - components successfully synthesized from an AOI
  sibling's averaged profile (see
  [`recover_missing_horizon_components()`](https://jjmaynard.github.io/soilSIM/reference/recover_missing_horizon_components.md)).

## Value

List with comprehensive metadata
