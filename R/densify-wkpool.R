## densify() for wkpool objects (hypertidy/wkpool >= 0.3.0)
##
## A wkpool is already bigcurve's native data model: a vertex pool plus
## directed segment index pairs, with segments as the atomic unit. The
## only translation needed is identity: wkpool references vertices by
## minted '.vx' ids (which survive subsetting and are not row positions),
## while the engine works in 0-based positions. So: ids -> positions on
## the way in, positions -> ids on the way out, with fresh ids minted for
## the vertices densification adds.
##
## Per-segment attributes (.feature) are carried through refinement via
## the engine's 'parent' index: every output segment inherits from the
## input segment it descends from. Added vertices are degree-2 by
## construction, so the arc-node structure of the pool is unchanged.

densify_wkpool <- function(x, target, source = "OGC:CRS84",
                           tolerance = NULL, pixels = 2048L,
                           max_depth = 16L) {
  if (!requireNamespace("wkpool", quietly = TRUE)) {
    stop("the 'wkpool' package is required to densify wkpool objects")
  }
  v <- wkpool::pool_vertices(x)
  if ("z" %in% names(v) || "m" %in% names(v)) {
    stop("densify() is strictly 2D: this pool carries z/m vertex values, ",
         "which cannot be carried through great circle densification. ",
         "Drop them (rebuild the pool from xy geometry) and try again.")
  }
  s <- wkpool::pool_segments(x)
  if (nrow(s) < 1L) {
    return(x)
  }

  ## minted .vx ids -> 0-based row positions for the engine
  p0 <- match(s$.vx0, v$.vx) - 1L
  p1 <- match(s$.vx1, v$.vx) - 1L
  stopifnot(!anyNA(p0), !anyNA(p1))

  tol <- tolerance %||%
    default_tolerance(cbind(v$x, v$y), target, source, pixels)

  out <- densify_mesh_cpp(v$x, v$y, p0, p1,
                          source, target, tol, as.integer(max_depth))

  ## mint ids for the appended vertices, continuing past the pool maximum
  n_old <- nrow(v)
  n_new <- length(out$x) - n_old
  vx_all <- c(v$.vx, if (n_new > 0L) max(v$.vx) + seq_len(n_new) else integer(0))
  vertices <- data.frame(.vx = vx_all, x = out$x, y = out$y)

  ## 0-based positions back to .vx ids; chain order preserves direction
  fields <- list(.vx0 = vx_all[out$s0 + 1L], .vx1 = vx_all[out$s1 + 1L])
  feat <- wkpool::pool_feature(x)
  if (!is.null(feat)) {
    fields$.feature <- as.integer(feat[out$parent + 1L])
  }

  ## wkpool does not (yet) export its constructor, so build the rcrd
  ## directly, honouring new_wkpool()'s invariants (checked above by
  ## construction: every .vx0/.vx1 is drawn from vertices$.vx)
  vctrs::new_rcrd(fields, pool = vertices, class = "wkpool")
}
