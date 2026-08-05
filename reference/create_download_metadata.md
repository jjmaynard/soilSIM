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
  validation_results
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

## Value

List with comprehensive metadata
