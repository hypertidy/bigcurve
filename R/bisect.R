#sf::sf_use_s2(FALSE)

#' @importFrom sf st_crs st_centroid st_transform st_linestring st_sfc st_distance st_point
#' @importFrom geosphere midPoint
mkdc <- function(x = "OGC:CRS84") {
  crs <- sf::st_crs("OGC:CRS84")
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
bisect <- function(x) {
  #sf::st_cast(lwgeom::st_split(x, sf::st_point(arc_centroid(x))))
  arc_split(x, arc_centroid(x))
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

