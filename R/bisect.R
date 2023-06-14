#' @importFrom terra project
rproj_xy <- function(x, target, ..., source = NULL) {
  if (is.null(source)) source <- "OGC:CRS84"
  l <- list(x = x[,1, drop = TRUE], y = x[,2, drop = TRUE])
  cpp_libproj_init_api()

  do.call(cbind, proj_coords(l, source, target))
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
  clnames <- function(x) {
    colnames(x) <- c("lon", "lat")
    x
  }
  ## here's the magic this is just projecting a point, comparing to its planar centroid
  geodist::geodist(clnames(rbind(arc_centroid(x),
                            rproj_xy( proj_centroid(x, proj), source = proj, target = default_crs()))), sequential = TRUE, measure = "vincenty")

}

arc_len <- function(x) {
  colnames(x)<- c("lon", "lat")
  geodist::geodist(x, sequential = TRUE, measure = "vincenty")
}

proj_centroid <- function(x, crs) {
  x <- rproj_xy(x, crs, source = default_crs())
  cbind(mean(x[,1L, drop = TRUE]),
        mean(x[,2L, drop = TRUE]))
}
## bigcurve originaly used geosphere but I think this is fine
## code is from SGAT
midpt <- function (p, fold = FALSE) {
  n <- nrow(p)
  rad <- pi/180
  p <- rad * p
  dlon <- diff(p[, 1L])
  lon1 <- p[-n, 1L]
  lat1 <- p[-n, 2L]
  lat2 <- p[-1L, 2L]
  bx <- cos(lat2) * cos(dlon)
  by <- cos(lat2) * sin(dlon)
  lat <- atan2(sin(lat1) + sin(lat2), sqrt((cos(lat1) + bx)^2 +
                                             by^2))/rad
  lon <- (lon1 + atan2(by, cos(lat1) + bx))/rad
#  if (fold)
 #   lon <- wrapLon(lon)
  cbind(lon, lat)
}

arc_centroid <- function(x){
  midpt(x)
}

arc_split <- function(x, pt) {
  rbind(x[1L, , drop = FALSE], pt, x[2L, , drop = FALSE])
}

bisect <- function(x, crs, dist = 1000, plt = FALSE) {
  cl <- curv_len(x, crs)
  al <- arc_len(x)
  if (!is.na(cl) && cl > 10000) {
    xx <- arc_split(x, arc_centroid(x))
    if (all(abs(x[1,] - xx[2, ]) < 0.001) || all(abs(x[2, ] - xx[2, ]) < 0.001 )) return(x)
    nn <<- nn + 1
    out <- rbind(bisect(xx[1:2L, , drop = FALSE], crs, dist, plt = plt),
                 bisect(xx[2:3L, , drop = FALSE], crs, dist, plt = plt))
    return(out)
  } else {
    x
  }
}

