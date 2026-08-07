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
