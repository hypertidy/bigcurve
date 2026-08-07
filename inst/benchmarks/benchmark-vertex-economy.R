# Vertex economy: adaptive densify() vs fixed-step segmentize
#
# Claim under test: for equal *measured* maximum deviation in the target
# projection, adaptive densification uses far fewer vertices than
# fixed-step densification, because fixed step must be sized for the
# worst local curvature and then spends that budget everywhere.
#
# All pipelines are judged by the same referee: sample each output
# segment's great circle finely, project, and measure the perpendicular
# distance to the segment's planar chord. This is bigcurve's own
# acceptance metric applied at high resolution, so it favours no one:
# it simply asks "how far does the drawn line stray from the true
# projected curve".

library(bigcurve)

## ---- referee: measured max deviation of a lon/lat path in target ----
max_deviation <- function(coords, target, source = "OGC:CRS84", k = 20) {
  p <- rproj_xy(coords, target, source = source)
  worst <- 0
  for (i in seq_len(nrow(coords) - 1L)) {
    ## k points along the great circle of this segment (slerp on sphere)
    a <- coords[i, ] * pi / 180; b <- coords[i + 1L, ] * pi / 180
    va <- c(cos(a[2]) * cos(a[1]), cos(a[2]) * sin(a[1]), sin(a[2]))
    vb <- c(cos(b[2]) * cos(b[1]), cos(b[2]) * sin(b[1]), sin(b[2]))
    ang <- acos(max(-1, min(1, sum(va * vb))))
    if (ang < 1e-12) next
    t <- seq(0, 1, length.out = k)
    m <- (sin((1 - t) * ang) %o% va + sin(t * ang) %o% vb) / sin(ang)
    ll <- cbind(atan2(m[, 2], m[, 1]), asin(pmax(-1, pmin(1, m[, 3])))) * 180 / pi
    q <- rproj_xy(ll, target, source = source)
    dx <- p[i + 1L, 1] - p[i, 1]; dy <- p[i + 1L, 2] - p[i, 2]
    d2 <- dx * dx + dy * dy
    if (!is.finite(d2) || d2 <= 0) next
    dz <- dx * (q[, 2] - p[i, 2]) - dy * (q[, 1] - p[i, 1])
    dev <- sqrt(dz * dz / d2)
    worst <- max(worst, dev[is.finite(dev)])
  }
  worst
}

count_vertices <- function(x) sum(vapply(x, nrow, integer(1)))

## ---- fixture: coarse world grid, the storage-minimal control mesh ----
p <- terra::as.polygons(terra::rast(terra::ext(-150, 150, -85, 85), res = 15))
g <- geos::as_geos_geometry(p)
pool <- wkpool::establish_topology(g)

## paths for the segmentize pipeline (sf wants geometry, not a pool)
sfx <- sf::st_as_sf(p)
sf::st_crs(sfx) <- "OGC:CRS84"

targets <- c(laea = "+proj=laea", stere = "+proj=stere +lat_0=-90",
             ortho = "+proj=ortho +lon_0=147 +lat_0=-42")
tol <- 5e4  # 50 km in projected metres

rows <- list()
for (nm in names(targets)) {
  tgt <- targets[[nm]]

  ## A. naive: no densification (the visibly broken baseline)
  naive <- lapply(seq_along(g), function(i) {
    xy <- wk::wk_coords(g[i])[, c("x", "y")]
    as.matrix(xy)
  })
  dev_naive <- max(vapply(naive, max_deviation, numeric(1), target = tgt))

  ## B. adaptive: densify the pool once at tol
  t_b <- system.time(
    dp <- densify(pool, tgt, tolerance = tol)
  )[["elapsed"]]
  nb <- nrow(wkpool::pool_vertices(dp))
  wb <- wkpool::segments_to_wkb(dp)
  dense_paths <- lapply(seq_along(wb), function(i)
    as.matrix(wk::wk_coords(wb[i])[, c("x", "y")]))
  dev_b <- max(vapply(dense_paths, max_deviation, numeric(1), target = tgt))

  ## C. fixed step: shrink dfMaxLength until measured deviation <= tol,
  ## i.e. give segmentize the step it truly needs for this projection
  step <- 2e6
  repeat {
    seg <- sf::st_segmentize(sfx, units::set_units(step, "m"))
    paths <- lapply(sf::st_geometry(seg), function(gm)
      do.call(rbind, lapply(unclass(gm), function(r) r[, 1:2])))
    dev_c <- max(vapply(paths, max_deviation, numeric(1), target = tgt))
    if (dev_c <= tol || step < 1e4) break
    step <- step / 2
  }
  nc <- count_vertices(paths)

  rows[[nm]] <- data.frame(
    target = nm, tol_m = tol,
    naive_dev_km = round(dev_naive / 1e3),
    adaptive_vertices = nb, adaptive_dev_km = round(dev_b / 1e3, 1),
    adaptive_secs = round(t_b, 3),
    fixed_step_km = step / 1e3, fixed_vertices = nc,
    fixed_dev_km = round(dev_c / 1e3, 1),
    vertex_ratio = round(nc / nb, 1)
  )
}
print(do.call(rbind, rows), row.names = FALSE)

## Optional, where applicable: s2's adaptive tessellator as honest prior
## art. It only speaks S2's built-in projections (e.g. plate carree,
## mercator), not arbitrary PROJ targets, which is the niche gap:
##   s2::s2_unprojection_filter / s2_projection_filter with
##   tessellate_tol, driven through wk::wk_handle.

## The other half of the story, once thin() exists: take an
## over-densified fixture (e.g. a 1-degree graticule, or the fixed-step
## output C above) and thin it back per arc at the same tolerance --
## the claim being that thin(C) lands close to B's vertex count, because
## the tolerance defines the budget and adaptive placement achieves it
## from either direction.
