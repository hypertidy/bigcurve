#sf::sf_use_s2(FALSE)

#' @importFrom sf st_crs st_centroid st_transform st_linestring st_sfc st_distance st_point
#' @importFrom geosphere midPoint
mkdc <- function(x = "OGC:CRS84") {
  crs <- sf::st_crs(x)
  function() {
    crs
  }
}
default_crs <- mkdc()

proj_centroid <- function(x, crs) {
  sf::st_centroid(sf::st_transform(x, crs))
}
arc_centroid <- function(x) {
  geosphere::midPoint(x[[1]][1, ], x[[1]][2, ])
}

arc_split <- function(x, pt) {
  sf::st_sfc(sf::st_linestring(rbind(x[[1]][1, ], pt)),
  sf::st_linestring(rbind(pt, x[[1]][2, ])), crs = sf::st_crs(x))
}
# bisect <- function(x) {
#   #sf::st_cast(lwgeom::st_split(x, sf::st_point(arc_centroid(x))))
#   arc_split(x, arc_centroid(x))
# }
#
bisect <- function(x, crs, dist = 1000) {
  cl <- curv_len(x, crs)
  if (cl > dist) {
    xx <- arc_split(x, arc_centroid(x))
    c(bisect(xx[1L], crs, dist), bisect(xx[2L], crs, dist))
  } else {
    x
  }
}
#' Segment
#'
#' @param x longitude (2, start point, end point)
#' @param y latitude (2, same as x)
#' @param crs angular projection to use, longlat defaults to 'OGC:CRS84'
#'
#' @return linestring sfc vector
#' @export
#' @importFrom stats runif
#' @examples
#' segment()
#' segment(c(0, 147), c(0, -42))
segment <- function(x = runif(2L, -180, 180), y = runif(2L, -90, 90), crs = default_crs()) {
  sf::st_sfc(sf::st_linestring(cbind(x, y)),crs = crs)
}

project_segment <- function(x, crs) {
  sf::st_transform(x, sf::st_crs(crs))
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

#' Curvature index for a given segment
#'
#' Calculate the perpendicular distance from the mid point of the projected
#' version of a segment to its arc mid-point.
#'
#' Yes this is the D3 Flawed Example logic.
#'
#' @param x linestring, single segment
#' @param proj projection string (proj, wkt, auth, whatever)
#'
#' @return distance in metres
#' @export
#'
#' @references
#' https://bost.ocks.org/mike/example/
#' @examples
#' s1 <- segment(c(0, 1), c(0, 1))
#' curv_len(s1, laea(147, -42))
#' curv_len(s1, laea(100, -30))
#' curv_len(s1, laea(70, -20))
#' curv_len(s1, laea(30, -10))
#' curv_len(s1, laea(0, 0))
curv_len <- function(x, proj) {
  ## here's the magic, unpick all this sf guff - this is just projecting a point, comparing to its planar centroid
  as.numeric(sf::st_distance(sf::st_transform(sf::st_sfc(sf::st_point(arc_centroid(x)), crs = default_crs()), proj),
                             proj_centroid(x, proj)))
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
mid_point <- function(x) {
  g <- geocentric(x)
  tg <- tail(g, -1)
  longlat(tg + (tg - head(g, -1))/2)
}

## for comparison, ultimately we should use a PROJ geocentric transform OR Karney's geodesic inverse code
# mid_point_proj <- function(x) {
#   # sf::sf_project("OGC:CRS84", "+proj=geocent +datum=WGS84",
#   #                cbind(x, 0))
#
#   p <- unclass(sf::st_transform(sf::st_sfc(sf::st_multipoint(cbind(x, 0)),
#                               crs = "OGC:CRS84"), "+proj=geocent +datum=WGS84")[[1]])
#   mp <- tail(p, -1) + (tail(p, -1) - head(p, -1))/2
#   unclass(sf::st_transform(sf::st_sfc(sf::st_multipoint(mp),
#                                       crs = "+proj=geocent +datum=WGS84"), "OGC:CRS84")[[1]])[,1:2, drop = FALSE]
# }

