# bigcurve 0.1.0

Complete rework of the package around a new adaptive resampling engine,
following d3-geo (Bostock). The prototype `bisect()` interface is gone.

## New user surface

* New exported function `densify()`, the single verb of the package.
  Accepts a two-column lon,lat matrix (a path), a list of such matrices,
  or a segment mesh (rgl mesh3d style `vb`/`is`), and returns the same
  kind of object with vertices inserted only where the target projection
  demands them.

* `rproj_xy()` is now exported: forward coordinate transformation with
  axis normalization (authority codes such as 'EPSG:4326' behave as
  lon,lat input) and `NA` rather than error for unprojectable points.

## Breaking changes

* Error is now controlled in the plane of the target projection: the
  `tolerance` argument is in projected units (usually metres), or
  derived from the projected extent as one part in `pixels` (default
  2048) when unset. The old `dist` metres-on-the-sphere argument (which
  was accepted but not honoured) is gone.

* Output coordinates are in the source coordinate system (geographic by
  default); project the result for plotting.

* Points outside the domain of the target projection never error:
  segments touching them are left un-densified and refinement stops at
  unprojectable midpoints. Clipping is out of scope.

* `bisect()`, `curv_len()`, `project_segment()` and `reproj_mesh3d()`
  are removed.

## Internals and dependencies

* New header-only C++ core with no R or PROJ dependency, driving
  refinement with a target-plane metric: perpendicular distance from
  the projected great circle midpoint to the planar chord, a parametric
  shear check, a forced split for arcs wider than 30 degrees, and a
  maximum depth. Each point is transformed exactly once and no inverse
  transforms are used, so refinement is dramatically faster than the
  0.0.x loop.

* Links directly to the system PROJ library (>= 8.0) via configure and
  pkg-config; the libproj package (removed from CRAN) is no longer
  used. Windows builds use the PROJ shipped with Rtools.

* Dropped terra, geodist, geosphere, glue, reproj, and all undeclared
  dependencies; hard dependencies are now cpp11 (LinkingTo) and system
  PROJ only.

* Spherical great circle midpoints are used (consistent with d3-geo);
  this is a display densification tool, not an ellipsoidal geodesy one.

* README rewritten around `densify()`, with a graticule example using
  `textures::segs()` (textures and maps added to Suggests).


# bigcurve 0.0.1

* First version with workable C++ components. Logic needs work but the pieces are in place. 


# bigcurve 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
