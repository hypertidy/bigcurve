test_that("bisection works", {
  s1 <- segment(c(0, 1), c(0, 1))
  expect_length(bisect(s1, "+proj=laea", 6378137), 4)
  expect_length(bisect(s1, "+proj=laea", 1), 8)
  expect_length(bisect(s1, "+proj=laea +lon_0=147", 1000), 8)

})

test_that("miscellaneous", {
  fun <- mkdc("EPSG:4267")
  expect_type(fun, "closure")


  expect_length(project_segment(segment(c(0, 1), c(0, 1)), "+proj=laea +lon_0=147"), 4)

  expect_equal(laea(), "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84")
})
