test_that("densify adds vertices where curvature demands them", {
  s <- segment(c(140, -164), c(-89, 80))
  crs <- "+proj=tmerc +lon_0=147 +lat_0=-42"
  d <- densify(s, crs)
  expect_true(is.matrix(d))
  expect_gt(nrow(d), 2L)
  ## endpoints preserved exactly
  expect_equal(unname(d[1L, ]), unname(s[1L, ]))
  expect_equal(unname(d[nrow(d), ]), unname(s[2L, ]))
  ## frame contract: this edge is written as crossing (140 -> -164,
  ## which unwraps to 196 in the start vertex's frame), so interior
  ## longitudes run continuously in [140, 196]; the final vertex keeps
  ## its written value (-164), the same sphere point
  interior <- d[-c(1L, nrow(d)), 1L]
  expect_true(all(interior >= 140 - 1e-9 & interior <= 196 + 1e-9))
})

test_that("straight-in-projection segments are left alone", {
  ## the equator is a straight line in an equator-centred laea
  s <- segment(c(-10, 10), c(0, 0))
  d <- densify(s, "+proj=laea +lon_0=0 +lat_0=0", tolerance = 1)
  expect_identical(nrow(d), 2L)
})

test_that("tolerance is monotone: tighter tolerance, more vertices", {
  s <- segment(c(140, -164), c(-89, 80))
  crs <- "+proj=laea +lon_0=147 +lat_0=-42"
  loose <- densify(s, crs, tolerance = 1e5)
  tight <- densify(s, crs, tolerance = 1e3)
  expect_gte(nrow(tight), nrow(loose))
  expect_gt(nrow(tight), 2L)
})

test_that("max_depth caps refinement", {
  s <- segment(c(140, -164), c(-89, 80))
  crs <- "+proj=laea +lon_0=147 +lat_0=-42"
  d1 <- densify(s, crs, tolerance = 1, max_depth = 2L)
  ## depth 2 per segment allows at most 3 interior vertices
  expect_lte(nrow(d1), 2L + 3L)
})

test_that("unprojectable regions do not error", {
  ## far side of an orthographic globe
  s <- segment(c(0, 10), c(0, 10))
  expect_silent(d <- densify(s, "+proj=ortho +lon_0=-170 +lat_0=0",
                             tolerance = 1000))
  expect_true(is.matrix(d))
})

test_that("mesh input round-trips with refined segments", {
  ## a tiny two-segment mesh
  mesh <- list(
    vb = rbind(c(140, -164, -89), c(-89, 80, 80), 0, 1),
    is = rbind(c(1L, 2L), c(2L, 3L))
  )
  crs <- "+proj=laea +lon_0=147 +lat_0=-42"
  out <- densify(mesh, crs, tolerance = 1e4)
  expect_true(ncol(out$vb) > ncol(mesh$vb))
  expect_true(ncol(out$is) > ncol(mesh$is))
  ## original vertices unchanged and still first
  expect_equal(out$vb[1:2, 1:3], mesh$vb[1:2, 1:3])
  ## every vertex referenced by 'is' exists
  expect_true(all(out$is >= 1L & out$is <= ncol(out$vb)))
})

test_that("list of matrices dispatches per element with shared tolerance", {
  s1 <- segment(c(140, -164), c(-89, 80))
  s2 <- segment(c(0, 60), c(-60, 10))
  crs <- "+proj=laea +lon_0=147 +lat_0=-42"
  out <- densify(list(s1, s2), crs)
  expect_length(out, 2L)
  expect_true(all(vapply(out, is.matrix, logical(1L))))
})

test_that("miscellaneous helpers", {
  expect_length(rproj_xy(segment(c(0, 1), c(0, 1)), "+proj=laea +lon_0=147"), 4L)
  expect_equal(laea(),
               "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84")
})

## a graticule as a segment mesh (rgl mesh3d style), shared vertices at
## every intersection: meridians every 'dx', parallels every 'dy'
graticule_mesh <- function(dx = 40, dy = 20, ext = c(-180, 180, -60, 60)) {
  lon <- seq(ext[1], ext[2], by = dx)
  lat <- seq(ext[3], ext[4], by = dy)
  g <- expand.grid(lon = lon, lat = lat)
  id <- matrix(seq_len(nrow(g)), length(lon))
  along <- cbind(c(id[-length(lon), ]), c(id[-1L, ]))
  up <- cbind(c(id[, -length(lat)]), c(id[, -1L]))
  structure(list(vb = rbind(g$lon, g$lat, 0, 1),
                 is = t(rbind(along, up))), class = "mesh3d")
}

test_that("a graticule mesh refines in one call", {
  prj <- "+proj=laea +lat_0=-60"
  mesh <- graticule_mesh()
  out <- densify(mesh, prj)
  expect_gt(ncol(out$vb), ncol(mesh$vb))
  expect_gt(ncol(out$is), ncol(mesh$is))
  ## input vertices preserved in place
  expect_equal(out$vb[1:2, seq_len(ncol(mesh$vb))], mesh$vb[1:2, ])
  ## refined index structure stays valid
  expect_true(all(out$is >= 1L & out$is <= ncol(out$vb)))
  ## and the result keeps its class
  expect_s3_class(out, class(mesh))
})

test_that("an edge repeated in a mesh is refined once and its vertices reused", {
  prj <- "+proj=laea +lon_0=50 +lat_0=30"
  mesh <- list(vb = rbind(c(0, 100), c(0, 50), 0, 1),
               is = cbind(c(1L, 2L), c(2L, 1L)))
  out <- densify(mesh, prj, tolerance = 2e4)
  n_new <- ncol(out$vb) - 2L
  expect_gt(n_new, 0L)
  ## two traversals, one set of minted vertices
  expect_identical(ncol(out$is), 2L * (n_new + 1L))
  fwd <- out$is[, seq_len(n_new + 1L)]
  bwd <- out$is[, n_new + 1L + seq_len(n_new + 1L)]
  expect_identical(bwd, fwd[2:1, rev(seq_len(ncol(fwd)))])
})

test_that("an edge refines to identical vertices in either direction", {
  crs <- "+proj=laea +lon_0=50 +lat_0=30"
  f <- densify(segment(c(0, 100), c(0, 50)), crs, tolerance = 2e4)
  r <- densify(segment(c(100, 0), c(50, 0)), crs, tolerance = 2e4)
  expect_gt(nrow(f), 2L)
  expect_identical(f, r[nrow(r):1, ])
})

test_that("parent maps every refined edge to its input edge", {
  ## engine-level contract used by the wkpool bridge
  x <- c(140, -164, -89)
  y <- c(-89, 80, 80)
  out <- densify_mesh_cpp(x, y, c(0L, 1L), c(1L, 2L),
                          "OGC:CRS84", "+proj=laea +lon_0=147 +lat_0=-42",
                          1e4, 16L)
  expect_length(out$parent, length(out$s0))
  expect_true(all(out$parent %in% c(0L, 1L)))
  ## chains are contiguous per parent, in input order
  expect_true(!is.unsorted(out$parent))
  ## each chain starts and ends at its input edge's endpoints
  expect_identical(out$s0[match(0L, out$parent)], 0L)
  expect_identical(out$s1[max(which(out$parent == 0L))], 1L)
  expect_identical(out$s0[match(1L, out$parent)], 1L)
  expect_identical(out$s1[max(which(out$parent == 1L))], 2L)
})

test_that("wkpool round trip: refine, preserve identity, carry features", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  ## two big adjacent polygons sharing the meridian at lon 80
  g <- wk::as_wkb(c(
    "POLYGON ((0 0, 80 0, 80 60, 0 60, 0 0))",
    "POLYGON ((80 0, 160 0, 160 60, 80 60, 80 0))"
  ))
  pool <- wkpool::establish_topology(g)
  crs <- "+proj=laea +lon_0=80 +lat_0=30"
  out <- densify(pool, crs, tolerance = 5e4)

  expect_s3_class(out, "wkpool")
  ## refinement happened
  expect_gt(length(out), length(pool))
  v0 <- wkpool::pool_vertices(pool)
  v1 <- wkpool::pool_vertices(out)
  expect_gt(nrow(v1), nrow(v0))
  ## original vertices first, ids and coordinates untouched
  expect_equal(v1[seq_len(nrow(v0)), ], v0)
  ## minted ids are new
  expect_true(all(v1$.vx[-seq_len(nrow(v0))] > max(v0$.vx)))
  ## every segment references a pooled vertex
  s1 <- wkpool::pool_segments(out)
  expect_true(all(s1$.vx0 %in% v1$.vx))
  expect_true(all(s1$.vx1 %in% v1$.vx))
  ## feature provenance carried, both features refined
  expect_setequal(unique(wkpool::pool_feature(out)),
                  unique(wkpool::pool_feature(pool)))
  ## node structure: rings traverse the shared meridian twice, and both
  ## traversals mint bitwise-identical vertices on it (all terms in the
  ## meridian midpoint are commutative sums), so merge_coincident fuses
  ## them into degree-4 vertices -- new nodes CAN appear, but only on
  ## the shared boundary at lon 80. The exact count is not asserted:
  ## forced splits leave children exactly on the 30-degree threshold,
  ## which is a strict float comparison and platform-sensitive.
  merged0 <- wkpool::merge_coincident(pool)
  merged1 <- wkpool::merge_coincident(out)
  n0 <- wkpool::find_nodes(merged0)
  n1 <- wkpool::find_nodes(merged1)
  expect_gte(length(n1), length(n0))
  vm <- wkpool::pool_vertices(merged1)
  new_nodes <- setdiff(n1, n0)
  if (length(new_nodes)) {
    expect_true(all(vm$x[match(new_nodes, vm$.vx)] == 80))
  }
})

test_that("added vertices are degree-2: exact on undupled topology", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  ## a linestring network has no duplicated traversals, so the node
  ## structure is exactly preserved: three lines meeting at a junction
  g <- wk::as_wkb(c(
    "LINESTRING (0 0, 60 40)",
    "LINESTRING (60 40, 120 0)",
    "LINESTRING (60 40, 60 85)"
  ))
  pool <- wkpool::merge_coincident(wkpool::establish_topology(g))
  out <- densify(pool, "+proj=laea +lon_0=60 +lat_0=40", tolerance = 5e4)
  expect_gt(nrow(wkpool::pool_vertices(out)), nrow(wkpool::pool_vertices(pool)))
  expect_identical(length(wkpool::find_nodes(out)),
                   length(wkpool::find_nodes(pool)))
})

test_that("wkpool with z vertices is refused, not mangled", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  g <- wk::as_wkb("LINESTRING Z (0 0 1, 80 60 2)")
  pool <- wkpool::establish_topology(g)
  skip_if(!"z" %in% names(wkpool::pool_vertices(pool)))
  expect_error(densify(pool, "+proj=laea +lon_0=40"), "strictly 2D")
})

test_that("antimeridian frame is preserved, east and west stay distinct", {
  crs <- "+proj=merc"
  ## an edge along lon 180 refines at 180
  e <- densify(segment(c(180, 180), c(-60, 60)), crs, tolerance = 5e4)
  expect_gt(nrow(e), 2L)
  expect_true(all(e[, 1L] == 180))
  ## its mirror along lon -180 refines at -180
  w <- densify(segment(c(-180, -180), c(-60, 60)), crs, tolerance = 5e4)
  expect_gt(nrow(w), 2L)
  expect_true(all(w[, 1L] == -180))
  ## and through the pool path: adjacent world halves keep both seams
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  g <- wk::as_wkb(c(
    "POLYGON ((0 -60, 180 -60, 180 60, 0 60, 0 -60))",
    "POLYGON ((-180 -60, 0 -60, 0 60, -180 60, -180 -60))"
  ))
  out <- densify(wkpool::establish_topology(g), crs, tolerance = 5e4)
  v <- wkpool::pool_vertices(out)
  expect_gt(sum(v$x == 180), 2L)   ## eastern seam has minted vertices
  expect_gt(sum(v$x == -180), 2L)  ## and so does the western seam
})

test_that("edges written as crossing get a locally continuous chain", {
  d <- densify(segment(c(170, -170), c(20, 20)), "+proj=merc",
               tolerance = 5e4)
  ## interior longitudes run monotonically through the frame [170, 190]
  interior <- d[-c(1L, nrow(d)), 1L]
  expect_true(all(diff(c(170, interior)) >= 0))
  expect_true(all(interior >= 170 & interior <= 190))
  ## bounded by the arc floor, not exploded to max_depth
  expect_lt(nrow(d), 200L)
})

test_that("arcs shorter than the tolerance are never refined", {
  ## ~22 km across the seam, tolerance 50 km: nothing to add
  d <- densify(segment(c(179.9, -179.9), c(0, 0)), "+proj=merc",
               tolerance = 5e4)
  expect_identical(nrow(d), 2L)
})

test_that("exactly-180-degree edges keep their written direction", {
  crs <- "+proj=merc"
  ## the four equator edges of a whole-world 2x2 grid, as rings write them
  e1 <- densify(segment(c(0, 180), c(0, 0)), crs, tolerance = 5e4)
  e2 <- densify(segment(c(180, 0), c(0, 0)), crs, tolerance = 5e4)
  w1 <- densify(segment(c(-180, 0), c(0, 0)), crs, tolerance = 5e4)
  w2 <- densify(segment(c(0, -180), c(0, 0)), crs, tolerance = 5e4)
  expect_true(all(e1[, 1L] >= 0 & e1[, 1L] <= 180))
  expect_true(all(e2[, 1L] >= 0 & e2[, 1L] <= 180))
  expect_true(all(w1[, 1L] >= -180 & w1[, 1L] <= 0))
  expect_true(all(w2[, 1L] >= -180 & w2[, 1L] <= 0))
  ## and they genuinely refined (forced splits on wide arcs)
  expect_gt(nrow(e1), 2L)
  expect_gt(nrow(w1), 2L)
})

test_that("the default tolerance works for a single parallel or meridian", {
  ## zero projected width in one direction
  d <- densify(segment(c(0, 100), c(60, 60)), "+proj=merc")
  expect_gt(nrow(d), 2L)
  m <- densify(segment(c(180, 180), c(-60, 60)), "+proj=laea +lon_0=150")
  expect_true(is.matrix(m))
})

test_that("a pool's own crs is used when 'source' is not given", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  g <- wk::wk_set_crs(wk::as_wkb("LINESTRING (0 0, 4000000 4000000)"),
                      "EPSG:3857")
  pool <- wkpool::establish_topology(g)
  out <- densify(pool, "+proj=laea +lat_0=45")
  v <- wkpool::pool_vertices(out)
  expect_gt(nrow(v), 2L)
  ## new vertices are in the pool's metres, between the endpoints
  new <- v[-(1:2), ]
  expect_true(all(new$x > 0 & new$x < 4e6))
  expect_identical(attr(out, "crs"), attr(pool, "crs"))
})

test_that("neighbouring polygons get identical vertices on a shared edge", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  ## the shared edge is a diagonal, not a meridian, so a direction-dependent
  ## midpoint would show up in the last bits
  g <- wk::as_wkb(c(
    "POLYGON ((0 0, 100 50, 0 60, 0 0))",
    "POLYGON ((0 0, 100 0, 100 50, 0 0))"
  ))
  pool <- wkpool::establish_topology(g)
  crs <- "+proj=laea +lon_0=50 +lat_0=30"
  out <- densify(pool, crs, tolerance = 2e4)
  v <- wkpool::pool_vertices(out)
  s <- wkpool::pool_segments(out)
  f <- wkpool::pool_feature(out)
  key <- function(ids) {
    vv <- v[match(unique(ids), v$.vx), ]
    vv <- vv[vv$.vx > max(wkpool::pool_vertices(pool)$.vx), ]
    paste(sprintf("%.17g", vv$x), sprintf("%.17g", vv$y))
  }
  k1 <- key(c(s$.vx0[f == 1L], s$.vx1[f == 1L]))
  k2 <- key(c(s$.vx0[f == 2L], s$.vx1[f == 2L]))
  ## vertices on the shared diagonal: lon/lat on the 0,0 -> 100,50 arc
  shared <- intersect(k1, k2)
  expect_gt(length(shared), 0L)
  ## merging fuses them: each minted vertex on the diagonal is shared
  m <- wkpool::merge_coincident(out)
  mv <- wkpool::pool_vertices(m)
  expect_lt(nrow(mv), nrow(v))
  ## and in a merged pool the shared edge is refined only once
  pm <- wkpool::merge_coincident(pool)
  om <- densify(pm, crs, tolerance = 2e4)
  n_minted <- nrow(wkpool::pool_vertices(om)) - nrow(wkpool::pool_vertices(pm))
  expect_identical(n_minted, nrow(v) - nrow(wkpool::pool_vertices(pool)) - length(shared))
})
