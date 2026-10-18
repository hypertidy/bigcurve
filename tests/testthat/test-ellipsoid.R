test_that("ellipsoid = TRUE puts new vertices on the WGS84 geodesic", {
  skip_if_not_installed("geographiclib")
  s <- rbind(c(147.3, -42.9), c(-0.1, 51.5))
  sph <- densify(s, "+proj=eqearth")
  ell <- densify(s, "+proj=eqearth", ellipsoid = TRUE)
  ## same vertex count, input vertices untouched
  expect_identical(nrow(ell), nrow(sph))
  expect_identical(ell[1L, ], s[1L, ])
  expect_identical(ell[nrow(ell), ], s[2L, ])
  ## on the geodesic: the two legs through each vertex add up to the whole
  s12 <- geographiclib::geodesic_inverse(s[1L, ], s[2L, ])$s12
  leg1 <- geographiclib::geodesic_inverse(s[rep(1L, nrow(ell)), ], ell)$s12
  leg2 <- geographiclib::geodesic_inverse(ell, s[rep(2L, nrow(ell)), ])$s12
  expect_lt(max(abs(leg1 + leg2 - s12)), 1e-3)
  ## and not on the great circle
  expect_gt(max(geographiclib::geodesic_distance(ell, sph)), 1e4)
})

test_that("ellipsoid replay is the same whichever way an edge is written", {
  skip_if_not_installed("geographiclib")
  s <- rbind(c(147.3, -42.9), c(-0.1, 51.5))
  fwd <- densify(s, "+proj=eqearth", ellipsoid = TRUE)
  rev <- densify(s[2:1, ], "+proj=eqearth", ellipsoid = TRUE)
  expect_identical(rev[nrow(rev):1, ], fwd)
})

test_that("ellipsoid replay keeps the longitude frame across 180", {
  skip_if_not_installed("geographiclib")
  x <- rbind(c(170, -40), c(-170, -45))
  d <- densify(x, "+proj=merc +lon_0=180", tolerance = 1000, ellipsoid = TRUE)
  inner <- d[-c(1L, nrow(d)), 1L]
  expect_true(all(inner > 170 & inner < 190))
})

test_that("ellipsoid replay keeps shared edges shared in a pool and a mesh", {
  skip_if_not_installed("geographiclib")
  skip_if_not_installed("wkpool")
  skip_if_not_installed("wk")
  g <- wk::as_wkb(c(
    "POLYGON ((0 0, 100 50, 0 60, 0 0))",
    "POLYGON ((0 0, 100 0, 100 50, 0 0))"
  ))
  pool <- wkpool::establish_topology(g)
  crs <- "+proj=laea +lon_0=50 +lat_0=30"
  out <- densify(pool, crs, tolerance = 2e4, ellipsoid = TRUE)
  pm <- wkpool::merge_coincident(pool)
  om <- densify(pm, crs, tolerance = 2e4, ellipsoid = TRUE)
  after <- wkpool::pool_vertices(wkpool::merge_coincident(out))
  first <- wkpool::pool_vertices(om)
  expect_identical(nrow(after), nrow(first))
  expect_identical(sort(complex(real = after$x, imaginary = after$y)),
                   sort(complex(real = first$x, imaginary = first$y)))

  mesh <- list(vb = rbind(c(0, 100), c(0, 50), 0, 1), is = matrix(c(1L, 2L), 2L))
  m <- densify(mesh, crs, tolerance = 2e4, ellipsoid = TRUE)
  p <- densify(rbind(c(0, 0), c(100, 50)), crs, tolerance = 2e4, ellipsoid = TRUE)
  expect_identical(sort(m$vb[1L, ]), sort(p[, 1L]))
})

test_that("ellipsoid must be TRUE or FALSE", {
  expect_error(densify(segment(c(0, 10), c(0, 10)), "+proj=merc", ellipsoid = NA))
})
