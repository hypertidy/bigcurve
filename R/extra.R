
mkdc <- function(x = "OGC:CRS84") {
  crs <- x
  function() {
    crs
  }
}
default_crs <- mkdc()

#' Segment
#'
#' @param x longitude (2, start point, end point)
#' @param y latitude (2, same as x)
#' @param crs angular projection to use, longlat defaults to 'OGC:CRS84'
#'
#' @return matrix
#' @export
#' @importFrom stats runif
#' @examples
#' segment()
#' segment(c(0, 147), c(0, -42))
segment <- function(x = runif(2L, -180, 180), y = runif(2L, -90, 90), crs = default_crs()) {
  cbind(x, y)
}

project_segment <- function(x, crs) {
  rproj_xy(x, crs)
}


#' Projection string
#'
#' Projection string, in old bad PROJ format.
#'
#' Note that easting/northing is very bad names for false X and Y offsets, but that's what they are called in PROJ.
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



geocentric <- function (x, rad = 6378137, exag = 1)
{
  cosLat = cos(x[, 2L] * pi/180)
  sinLat = sin(x[, 2L] * pi/180)
  cosLon = cos(x[, 1L] * pi/180)
  sinLon = sin(x[, 1L] * pi/180)
  cbind(rad * cosLat * cosLon,
        rad * cosLat * sinLon,
        rad * sinLat)
}
clamp1 <- function(x) {
  xgt <- x > 1
  xlt <- x < -1
  if (any(xgt)) x[xgt] <- 1
  if (any(xlt)) x[xlt] <- -1
  x
}
longlat <- function(x, rad = 6378137) {
  Zr <- clamp1(x[,3]/rad)
  cbind(atan2 (x[,2L], x[,1L]), asin (Zr)) * 180/pi
}
#' @importFrom utils head tail
mid_point <- function(x) {
  g <- geocentric(x)
  longlat(cbind(mean(g[,1]), mean(g[,2]), mean(g[,3])))
}



