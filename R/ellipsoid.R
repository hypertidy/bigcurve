## Ellipsoidal replay of a spherical densification (ellipsoid = TRUE)
##
## The engine splits every arc at its great circle midpoint, so each new
## vertex sits at a dyadic fraction t (k / 2^d) of its input edge's arc.
## The sphere decides how many vertices an edge needs and where along it
## they go; this moves each new vertex to the same fraction of the WGS84
## geodesic between the same two endpoints, using the geographiclib
## package. Input vertices are never moved.
##
## Identity: the fraction is measured from the edge's canonical start (the
## endpoint with the smaller lon, lat, as in the engine) and snapped to the
## dyadic grid, so an edge refines to the same vertices bit for bit
## whichever way it is written.

check_ellipsoid <- function() {
  if (!requireNamespace("geographiclib", quietly = TRUE)) {
    stop("'ellipsoid = TRUE' needs the 'geographiclib' package")
  }
}

## angle between lon/lat rows on the unit sphere, in radians
arc_angle <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  x1 <- cos(lat1 * r) * cos(lon1 * r); y1 <- cos(lat1 * r) * sin(lon1 * r)
  z1 <- sin(lat1 * r)
  x2 <- cos(lat2 * r) * cos(lon2 * r); y2 <- cos(lat2 * r) * sin(lon2 * r)
  z2 <- sin(lat2 * r)
  cx <- y1 * z2 - z1 * y2; cy <- z1 * x2 - x1 * z2; cz <- x1 * y2 - y1 * x2
  atan2(sqrt(cx * cx + cy * cy + cz * cz), x1 * x2 + y1 * y2 + z1 * z2)
}

## out: densify_mesh_cpp() result; n_old: number of input vertices;
## s0, s1: 0-based input edges. Returns out with new vertices moved.
ellipsoid_replay <- function(out, n_old, s0, s1, max_depth) {
  n_new <- length(out$x) - n_old
  if (n_new < 1L) return(out)
  new_id <- n_old + seq_len(n_new) - 1L

  ## parent input edge of each new vertex (each one starts an output segment)
  k <- match(new_id, out$s0)
  p <- out$parent[k] + 1L
  a <- s0[p] + 1L
  b <- s1[p] + 1L
  ax <- out$x[a]; ay <- out$y[a]; bx <- out$x[b]; by <- out$y[b]
  rev <- bx < ax | (bx == ax & by < ay)
  px <- ifelse(rev, bx, ax); py <- ifelse(rev, by, ay)
  qx <- ifelse(rev, ax, bx); qy <- ifelse(rev, ay, by)

  vx <- out$x[new_id + 1L]; vy <- out$y[new_id + 1L]
  grid <- 2^min(max_depth, 50L)
  t <- round(arc_angle(px, py, vx, vy) / arc_angle(px, py, qx, qy) * grid) / grid

  ## one inverse problem per distinct canonical edge
  key <- paste(px, py, qx, qy)
  u <- !duplicated(key)
  inv <- geographiclib::geodesic_inverse(cbind(px[u], py[u]), cbind(qx[u], qy[u]))
  j <- match(key, key[u])
  d <- geographiclib::geodesic_direct(cbind(px, py), inv$azi1[j], t * inv$s12[j])

  ## keep the longitude frame the engine gave the vertex
  lon <- d$lon2 + 360 * round((vx - d$lon2) / 360)
  out$x[new_id + 1L] <- lon
  out$y[new_id + 1L] <- d$lat2
  out
}

## a path through the mesh engine, so each new vertex knows its edge
densify_path_ellipsoid <- function(x, source, target, tol, max_depth) {
  n <- nrow(x)
  if (n < 2L) return(x[, 1:2, drop = FALSE])
  s0 <- seq_len(n - 1L) - 1L
  s1 <- s0 + 1L
  out <- densify_mesh_cpp(x[, 1L, drop = TRUE], x[, 2L, drop = TRUE],
                          s0, s1, source, target, tol, max_depth)
  out <- ellipsoid_replay(out, n, s0, s1, max_depth)
  idx <- c(out$s0, out$s1[length(out$s1)]) + 1L
  cbind(out$x[idx], out$y[idx])
}
