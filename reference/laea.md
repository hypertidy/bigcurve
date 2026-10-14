# Projection string

Lambert azimuthal equal area projection string, in old-style PROJ
format.

## Usage

``` r
laea(lon_0 = 0, lat_0 = 0, x_0 = 0, y_0 = 0)
```

## Arguments

- lon_0:

  centre longitude

- lat_0:

  centre latitude

- x_0:

  false easting (default 0 is fine)

- y_0:

  false northing (default 0 is fine)

## Value

character string, projection

## Details

Note that easting/northing are poor names for false X and Y offsets, but
that is what they are called in PROJ.

## Examples

``` r
laea()
#> [1] "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84"
laea(147, -42)
#> [1] "+proj=laea +lon_0=147 +lat_0=-42 +x_0=0 +y_0=0 +datum=WGS84"
```
