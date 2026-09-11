# Tag a list of percentile-value rasters for `fuse_property()`'s raster method

[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md)
dispatches on the class of `prior`: a plain numeric vector or plain list
goes to
[`fuse_property.default()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md)
(scalar/vector fusion); a list tagged with this class goes to
[`fuse_property.percentile_rasters()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md)
(raster-native fusion, forwarding to
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
unchanged). A plain, unclassed
[`list()`](https://rdrr.io/r/base/list.html) of percentile-value rasters
cannot be distinguished from a closed-form parameter list by class
alone, hence the explicit tag.

## Usage

``` r
new_percentile_rasters(x)
```

## Arguments

- x:

  A list of
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  percentile-value layers.

## Value

`x`, with class `c("percentile_rasters", class(x))`.
