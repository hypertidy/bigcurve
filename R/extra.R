mkdc <- function(x = "EPSG:4326") {
  crs <- x
  function() {
    crs
  }
}
default_crs <- mkdc()

#' Segment
#'
#' A two-point path as a two-column matrix, by default a random one.
#'
#' @param x longitude (2, start point, end point)
#' @param y latitude (2, same as x)
#'
#' @return two-column matrix
#' @export
#' @importFrom stats runif
#' @examples
#' segment()
#' segment(c(0, 147), c(0, -42))
segment <- function(x = runif(2L, -180, 180), y = runif(2L, -90, 90)) {
  cbind(x, y)
}

#' Projection string
#'
#' Lambert azimuthal equal area projection string, in old-style PROJ format.
#'
#' Note that easting/northing are poor names for false X and Y offsets, but
#' that is what they are called in PROJ.
#' @param lon_0 centre longitude
#' @param lat_0 centre latitude
#' @param x_0 false easting (default 0 is fine)
#' @param y_0 false northing (default 0 is fine)
#'
#' @return character string, projection
#' @export
#'
#' @examples
#' laea()
#' laea(147, -42)
laea <- function(lon_0 = 0, lat_0 = 0, x_0 = 0, y_0 = 0) {
  glue::glue("+proj=laea +lon_0={lon_0} +lat_0={lat_0} +x_0={x_0} +y_0={y_0} +datum=WGS84")
}
