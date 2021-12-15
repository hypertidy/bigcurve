
path_to_linestrings <-function(x) {
  sc <- silicate::SC0(x)
  #g <- silicate::sc_path(x)
  v <- silicate::sc_vertex(sc)
  idx <- do.call(rbind, lapply(sc$object$topology_, \(d) as.matrix(d[c(".vx0", ".vx1")])))

  v <- v[c(t(idx)), ]
  v$segment_id <- rep(seq_len(dim(idx)[1L]), each = 2L)
  sf::st_set_crs(sfheaders::sfc_linestring(v, linestring_id = "segment_id"), sf::st_crs(x))
}

sf::sf_use_s2(FALSE)
x <- sf::st_crop(subset(rnaturalearth::ne_countries(returnclass = "sf"), sovereignt == "Antarctica"),
                 sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90)))




g <- silicate::sc_path(x)
l <- path_to_linestrings(x)

plot(l)
proj <- "+proj=vandg"
plot(sf::st_transform(l, proj))



l1 <- vector("list", length(l))
## these are in metres, via lwgeom
lens <- as.numeric(sf::st_length(l))
ex <-  sf::st_bbox(sf::st_transform(l, proj))[c(1, 3, 2, 4)]
minr <- min(c(diff(ex[1:2]), diff(ex[3:4])))
dist <- 50000
for (i in seq_along(l)) {
#  rat <- minr/lens
  if (lens[i] > 5e5) {
    l1[[i]] <- bisect(l[i], proj, dist)
    print(i)
  } else {
    l1[[i]] <- l[i]
  }
}


plot(sf::st_transform(sf::st_sfc(unlist(l1, recursive = FALSE), crs = sf::st_crs(x)), proj))



