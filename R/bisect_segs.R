bisect_segs_cpp1 <- function(segmesh, prj, src = "OGC:CRS84", threshold = NULL) {
  nlast <- 0
  while(TRUE) {
    d <- bisect_cpp(list(xxs(segmesh), yys(segmesh)), src, prj)

    if (is.null(threshold)) {
      bad <- is.na(d$x) | is.na(d$y)
      bb <- apply(terra::project(cbind(d$x, d$y)[!bad, ], to = prj, from = src), 2, range)
      threshold <- min(apply(bb, 2, diff) / 4000)
    }
    idx <- which(d$dist > threshold)
    if (length(idx) == nlast) break; ## avoid infinite loop
    if (length(idx) < 1) break;
    #verts <- segmesh$is[,idx]
    nv <- ncol(segmesh$vb)
    ## add in all the new vertices

    for (i in seq_along(idx)) {
      if (nlast == 0 && i == 1) print(c(d$x[idx[i]], d$y[idx[i]], 0, 1))
      segmesh$vb <- cbind(segmesh$vb, c(d$x[idx[i]], d$y[idx[i]], 0, 1))
      segmesh$is <- cbind(segmesh$is, c(segmesh$is[1,idx[i]], nv + i), ## first newly bisected segmented
                          c(nv + i, segmesh$is[2,idx[i]])  ## second ""
      )
    }
    ## now that we've added new vertices and segments, we can safely delete the original segments identified for bisection
    segmesh$is <- segmesh$is[,-idx]
    nlast <- length(idx)
  }
  segmesh
}

xxs <- function(x) {
  x$vb[1L, x$is]
}
yys <- function(x) {
  x$vb[2L, x$is]
}

bisect_segs_r1 <- function(segmesh, prj, src = "OGC:CRS84", threshold = NULL) {



  nlast <- 0
  while(TRUE) {
    mp_gc <- do.call(cbind, mid_pt_pairs_gc(xxs(segmesh), yys(segmesh)))
    xy <- terra::project(cbind(xxs(segmesh), yys(segmesh)), to = prj, from = src)
    mp <- terra::project(do.call(cbind, mid_pt_pairs_plane(xy[,1], xy[,2])), from = prj, to = src)
    colnames(mp_gc) <- c("lon", "lat")
    colnames(mp) <- c("lon", "lat")

    d <- geosphere::distGeo(mp_gc, mp)

    if (is.null(threshold)) {
      bad <- is.na(xy[,1]) | is.na(xy[,2])

      bb <- apply(xy[!bad, ], 2, range)
      threshold <- min(apply(bb, 2, diff) / 4000, na.rm = TRUE)
    }
    idx <- which(d > threshold)
    if (nlast > 5000) break;
    if (length(idx) == nlast) break; ## avoid infinite loop
    if (length(idx) < 1) break;
    #verts <- segmesh$is[,idx]
    nv <- ncol(segmesh$vb)
    ## add in all the new vertices
#browser()
    for (i in seq_along(idx)) {
      if (nlast == 0 && i == 1) print(c(mp_gc[idx[i], ], 0, 1))
      segmesh$vb <- cbind(segmesh$vb, c(mp_gc[idx[i], ], 0, 1))
      ## first newly bisected segmented
      segmesh$is <- cbind(segmesh$is, c(segmesh$is[1,idx[i]], nv + i),
                          c(nv + i, segmesh$is[2,idx[i]])  ## second ""
      )
    }
    plot(segmesh)
    scan("", 1)
    ## now that we've added new vertices and segments, we can safely delete the original segments identified for bisection
    #segmesh$is <- segmesh$is[,-idx]
    nlast <- length(idx)
    print(i)
  }
  segmesh
}
