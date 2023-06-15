mid_point_proj <- function(x, proj) {
  out <- rproj_xy(x, proj)
  cbind(mean(out[,1]), mean(out[,2]))
}

prs <- function (x)
{
  cbind(head(x, -1), tail(x, -1))
}

split_n <- function(x, n = 1) {
  split(x, rep(seq_len(length(x)/2), each = n))
}

SC0_c <- function(x, y, crs = NA_character_) {
  verts <- tibble::tibble(x_ = x, y_ = y)

  uj <- unjoin::unjoin(verts, "x_", "y_")
  splits <- split_n(uj$data$.idx0, n = 2)
  lt <- lapply(seq_along(splits), function(.x) tibble::tibble(.vx0 = splits[[.x]][1], .vx1 = splits[[.x]][2], path_ = .x))

  structure(list(object = tibble::tibble(object_ = 1, topology_ = lt),
                 vertex = uj$.idx0[order(uj$.idx0$.idx0), c("x_", "y_")], meta = tibble::tibble(proj = crs, ctime = Sys.time())), class = c("SC0", "sc"))
}


mk_mesh3d <- function(x, y) {
  uj <- unjoin::unjoin(tibble::tibble(x = x, y = y), "x", "y")
  rgl::mesh3d(vertices = rbind(uj$.idx0$x, uj$.idx$y, 0, 1)[,order(uj$.idx$.idx0),  drop = FALSE],
              segments = matrix(uj$data$.idx0, nrow = 2L))
}


#' @importFrom reproj reproj reproj_xy
#' @export
reproj_mesh3d <- function(x, target, ..., source = NULL) {
  x$vb[1:2, ] <- t(reproj_xy(t(x$vb[1:2, ]), target, source = source))
  x
}
