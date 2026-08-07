#' Transform coordinates
#'
#' Forward-transform a two-column matrix of coordinates with PROJ. Authority
#' codes are axis-normalized, so input is always x,y (lon,lat) order and
#' output is easting,northing - even for 'EPSG:4326'.
#'
#' Points outside the domain of `target` are returned as `NA` rather than
#' raising an error.
#'
#' @param x two-column matrix of coordinates in `source`
#' @param target output coordinate system (anything PROJ accepts)
#' @param ... ignored
#' @param source input coordinate system, default 'OGC:CRS84'
#'
#' @return two-column matrix of transformed coordinates
#' @export
#' @examples
#' rproj_xy(cbind(147, -42), "+proj=laea +lon_0=147 +lat_0=-42")
rproj_xy <- function(x, target, ..., source = NULL) {
  if (is.null(source)) source <- "OGC:CRS84"
  l <- list(x = x[, 1L, drop = TRUE], y = x[, 2L, drop = TRUE])
  out <- do.call(cbind, proj_coords(l, source, target))
  out[!is.finite(out)] <- NA_real_
  out
}
