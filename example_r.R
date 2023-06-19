library(textures)
library(bigcurve)

EXT <- c(-180, 180,  -65, 65)

xxs <- function(x) {
  x$vb[1L, x$is]
}
yys <- function(x) {
  x$vb[2L, x$is]
}

prj <- "+proj=lcc +lat_1=-60 +lat_2=-30 +lon_0=100"
#prj <- "+proj=ortho +lon_0=100"
#prj <- "+proj=stere +lat_0=-90 +lon_0=147 +lat_ts=-71"
#prj <- "+proj=tmerc +R=6378137"
prj <- "+proj=laea +lat_0=-60"
#prj <- "+proj=ortho +lon_0=-35.000000 +lat_0=58.298038 +R=6378137"
#prj <- "+proj=tobmerc"
src <- "+proj=longlat +R=6378137"
bigcurve:::cpp_libproj_init_api()

segmesh <- segs(c(40, 20), extent = EXT)


#range(d)
#segmesh <- segmesh2

xy <- rproj_xy(cbind(xxs(segmesh), yys(segmesh)), prj)


nlast <- 0
#ran <- function(x) {x <- x[!is.na(x)]; x <- x[is.finite(x)]; range(x)}
while(TRUE) {

  mp_gc <- do.call(cbind, mid_pt_pairs_gc(xxs(segmesh), yys(segmesh)))
  xy <- rproj_xy(cbind(xxs(segmesh), yys(segmesh)), prj)
  mp <- rproj_xy(do.call(cbind, mid_pt_pairs_plane(xy[,1], xy[,2])), "OGC:CRS84", source = prj)
  d <- list(x = mp_gc[,1], y = mp_gc[,2], dist = geosphere::distGeo(mp_gc, mp))


#bb <- apply(terra::project(xy, to = src, from = prj), 2, range, na.rm = TRUE)

  bb <- apply(xy, 2, range, na.rm = TRUE)

idx <- which(d$dist > (min(apply(bb, 2, diff) / 2048)))
if (length(idx) == nlast) break; ## avoid infinite loop
if (length(idx) < 1) break;
#verts <- segmesh$is[,idx]
nv <- ncol(segmesh$vb)
## add in all the new vertices

for (i in seq_along(idx)) {
  segmesh$vb <- cbind(segmesh$vb, c(d$x[idx[i]], d$y[idx[i]], 0, 1))
  segmesh$is <- cbind(segmesh$is, c(segmesh$is[1,idx[i]], nv + i), ## first newly bisected segmented
                                  c(nv + i, segmesh$is[2,idx[i]])  ## second ""
  )
}
## now that we've added new vertices and segments, we can safely delete the original segments identified for bisection
segmesh$is <- segmesh$is[,-idx]

print(length(idx))
nlast <- length(idx)
}

projsegmesh <- segmesh
projsegmesh$vb[1:2, ] <- t(terra::project(t(segmesh$vb[1:2, ]), to = prj, from = src))
plot(projsegmesh)
points(t(projsegmesh$vb[1:2, ]), pch = 19, cex = 0.5, col = "red")

## segment length
x0 <- t(projsegmesh$vb[1:2, projsegmesh$is[1, ]])
x1 <- t(projsegmesh$vb[1:2, projsegmesh$is[2, ]])

dx <- x0[,1] - x1[,1]
dy <- x0[,2] - x1[,2]
len <- sqrt(dx * dx + dy * dy)
#bad <- len > 1e7
#projsegmesh$is <- projsegmesh$is[,!bad]
plot(projsegmesh, asp = 1)
points(t(projsegmesh$vb), pch = 19, cex = .4)



xy <- wk::xy(projsegmesh$vb[1, ], projsegmesh$vb[2, ])
x <- geos::geos_polygonize(wk::wk_collection(do.call(c, lapply(split(segmesh$is, rep(1:ncol(segmesh$is), each = 2)), function(.x) wk::wk_linestring(xy[.x])))))
x <- wk::wk_flatten(x)
plot(x, col = hcl.colors(length(x)))
