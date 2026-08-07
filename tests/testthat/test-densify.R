test_that("densify adds vertices where curvature demands them", {
  s <- segment(c(140, -164), c(-89, 80))
  crs <- "+proj=tmerc +lon_0=147 +lat_0=-42"
  d <- densify(s, crs)
  expect_true(is.matrix(d))
  expect_gt(nrow(d), 2L)
  ## endpoints preserved exactly
  expect_equal(unname(d[1L, ]), unname(s[1L, ]))
  expect_equal(unname(d[nrow(d), ]), unname(s[2L, ]))
  ## longitudes wrapped
  expect_true(all(d[, 1L] >= -180 & d[, 1L] < 180))
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

test_that("a textures::segs graticule refines in one call", {
  skip_if_not_installed("textures")
  prj <- "+proj=laea +lat_0=-60"
  mesh <- textures::segs(c(40, 20), extent = c(-180, 180, -65, 65))
  out <- densify(mesh, prj)
  expect_gt(ncol(out$vb), ncol(mesh$vb))
  expect_gt(ncol(out$is), ncol(mesh$is))
  ## input vertices preserved in place
  expect_equal(out$vb[1:2, seq_len(ncol(mesh$vb))], mesh$vb[1:2, ])
  ## refined index structure stays valid
  expect_true(all(out$is >= 1L & out$is <= ncol(out$vb)))
  ## and the result is still plottable as the same class
  expect_s3_class(out, class(mesh))
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
  ## added vertices are degree-2: node structure is unchanged
  merged0 <- wkpool::merge_coincident(pool)
  merged1 <- wkpool::merge_coincident(out)
  expect_identical(length(wkpool::find_nodes(merged1)),
                   length(wkpool::find_nodes(merged0)))
})

test_that("wkpool with z vertices is refused, not mangled", {
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  g <- wk::as_wkb("LINESTRING Z (0 0 1, 80 60 2)")
  pool <- wkpool::establish_topology(g)
  skip_if(!"z" %in% names(wkpool::pool_vertices(pool)))
  expect_error(densify(pool, "+proj=laea +lon_0=40"), "strictly 2D")
})
