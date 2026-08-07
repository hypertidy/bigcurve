`%||%` <- function(a, b) if (is.null(a)) b else a

#' Adaptively densify lines for a projection
#'
#' Insert vertices along great circle arcs exactly where a target projection
#' bends them, so that straight chords drawn between output vertices stay
#' within `tolerance` of the true projected curve. This is the adaptive
#' resampling scheme of d3-geo (Bostock), generalized to paths and to
#' segment meshes.
#'
#' The error metric is evaluated in the plane of the `target` projection:
#' the great circle midpoint of each candidate segment is forward-projected
#' and compared to the segment's planar chord (perpendicular distance, plus
#' a parametric check that catches strong shear). Segments are split
#' recursively until they pass or `max_depth` is reached. No inverse
#' transforms are used and each point is projected exactly once.
#'
#' Coordinates are returned in the coordinates of `source` (geographic by
#' default): project the result yourself for plotting, e.g. with
#' [rproj_xy()] or your tool of choice.
#'
#' Points that cannot be projected (for example beyond the horizon of an
#' orthographic projection) do not error: segments touching them are left
#' un-densified, and refinement stops at unprojectable midpoints. Clipping
#' is out of scope here.
#'
#' @param x a two-column matrix of lon,lat (a path), a list of such
#'   matrices, or a segment mesh: a list with `vb` (4 x n vertex matrix,
#'   rows x, y, z, h) and `is` (2 x m segment index matrix, 1-based), as
#'   used by rgl's mesh3d
#' @param target projection for which to densify (proj string, WKT,
#'   authority code - anything PROJ accepts)
#' @param source coordinate system of the input, default 'OGC:CRS84'
#'   (longitude, latitude); authority codes are axis-normalized so
#'   'EPSG:4326' input is still given as lon,lat
#' @param tolerance maximum allowed deviation, in units of `target`
#'   (usually metres); if `NULL` (the default) it is derived from the
#'   projected extent of the input as `min(diff(range))/pixels`
#' @param pixels notional output width/height used to derive the default
#'   `tolerance`, i.e. "keep error under about one part in `pixels` of the
#'   plot"; ignored when `tolerance` is supplied
#' @param max_depth maximum bisection depth per input segment (each level
#'   at most doubles the vertex count of a segment)
#'
#' @return the same kind of object as `x`, with vertices added: a matrix
#'   for a path, a list of matrices for a list, a mesh with appended
#'   vertices and refined `is` for a mesh
#' @export
#'
#' @references <https://bost.ocks.org/mike/example/>
#' @examples
#' s <- segment(c(140, -164), c(-89, 80))
#' d <- densify(s, "+proj=laea +lon_0=147 +lat_0=-42")
#' nrow(d)
#' plot(rproj_xy(d, "+proj=laea +lon_0=147 +lat_0=-42"), type = "l", asp = 1)
densify <- function(x, target, source = "OGC:CRS84",
                    tolerance = NULL, pixels = 2048L, max_depth = 16L) {
  stopifnot(is.numeric(pixels), pixels >= 1, is.numeric(max_depth))
  max_depth <- as.integer(max_depth)

  if (is.matrix(x)) {
    tol <- tolerance %||% default_tolerance(x, target, source, pixels)
    out <- densify_path_cpp(x[, 1L, drop = TRUE], x[, 2L, drop = TRUE],
                            source, target, tol, max_depth)
    return(cbind(out$x, out$y))
  }

  if (is.list(x) && !is.null(x$vb) && !is.null(x$is)) {
    verts <- t(x$vb[1:2, , drop = FALSE])
    tol <- tolerance %||% default_tolerance(verts, target, source, pixels)
    out <- densify_mesh_cpp(verts[, 1L], verts[, 2L],
                            as.integer(x$is[1L, ]) - 1L,
                            as.integer(x$is[2L, ]) - 1L,
                            source, target, tol, max_depth)
    x$vb <- rbind(out$x, out$y, 0, 1)
    x$is <- rbind(out$s0 + 1L, out$s1 + 1L)
    return(x)
  }

  if (is.list(x) && all(vapply(x, is.matrix, logical(1L)))) {
    tol <- tolerance %||%
      default_tolerance(do.call(rbind, x), target, source, pixels)
    return(lapply(x, densify, target = target, source = source,
                  tolerance = tol, pixels = pixels, max_depth = max_depth))
  }

  stop("'x' must be a matrix, a list of matrices, or a segment mesh (vb/is)")
}

## default tolerance: one part in 'pixels' of the smaller side of the
## projected bounding box of the input vertices (the d3 screen-space idea)
default_tolerance <- function(verts, target, source, pixels) {
  xy <- rproj_xy(verts, target, source = source)
  ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  if (!any(ok)) {
    stop("no input vertex is projectable in 'target'; supply 'tolerance' explicitly")
  }
  side <- min(diff(range(xy[ok, 1L])), diff(range(xy[ok, 2L])))
  if (!is.finite(side) || side <= 0) {
    stop("degenerate projected extent; supply 'tolerance' explicitly")
  }
  side / pixels
}
