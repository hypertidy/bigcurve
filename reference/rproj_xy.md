# Transform coordinates

Forward-transform a two-column matrix of coordinates with PROJ.
Authority codes are axis-normalized, so input is always x,y (lon,lat)
order and output is easting,northing - even for 'EPSG:4326'.

## Usage

``` r
rproj_xy(x, target, ..., source = NULL)
```

## Arguments

- x:

  two-column matrix of coordinates in `source`

- target:

  output coordinate system (anything PROJ accepts)

- ...:

  ignored

- source:

  input coordinate system, default 'EPSG:4326' (axis-normalized, so
  input is lon,lat regardless of authority order)

## Value

two-column matrix of transformed coordinates

## Details

Points outside the domain of `target` are returned as `NA` rather than
raising an error.

## Examples

``` r
rproj_xy(cbind(147, -42), "+proj=laea +lon_0=147 +lat_0=-42")
#>      x y
#> [1,] 0 0
```
