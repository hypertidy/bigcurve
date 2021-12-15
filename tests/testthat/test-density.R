test_that("bisection works", {
  s1 <- segment(c(0, 1), c(0, 1))
  expect_length(bisect(s1, "+proj=laea", 6378137), 1)
  expect_length(bisect(s1, "+proj=laea", 1), 2)
  expect_length(bisect(s1, "+proj=laea +lon_0=147", 1000), 2)

})

test_that("miscellaneous", {
  fun <- mkdc("EPSG:4267")
  expect_type(fun, "closure")
  expect_equal(fun(), sf::st_crs(4267))

  expect_s3_class(project_segment(segment(c(0, 1), c(0, 1)), "+proj=laea +lon_0=147"), "sfc")

  expect_equal(laea(), "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84")
})
