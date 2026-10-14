# Segment

A two-point path as a two-column matrix, by default a random one.

## Usage

``` r
segment(x = runif(2L, -180, 180), y = runif(2L, -90, 90))
```

## Arguments

- x:

  longitude (2, start point, end point)

- y:

  latitude (2, same as x)

## Value

two-column matrix

## Examples

``` r
segment()
#>              x         y
#> [1,] -150.9300  18.13696
#> [2,]  120.3599 -61.70248
segment(c(0, 147), c(0, -42))
#>        x   y
#> [1,]   0   0
#> [2,] 147 -42
```
