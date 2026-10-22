# bigcurve: Adaptive Densification of Lines for Map Projections

Densify lines adaptively, as required by the curvature a map projection
imposes on them. Edges are read as great circle arcs, and vertices are
inserted along them only where the projected curve deviates from a
straight chord by more than a tolerance measured on the map, following
the adaptive resampling method of the 'd3-geo' 'JavaScript' library
described in Bostock (2016) <https://github.com/d3/d3-geo>. Projections
are computed with the 'PROJ' coordinate transformation library
(<https://proj.org/>), so any coordinate reference system it supports
can be a target. Inputs can be coordinate matrices, segment meshes or
vertex pool topology from 'wkpool', and edges shared by neighbouring
polygons receive identical new vertices.

## See also

Useful links:

- <https://github.com/hypertidy/bigcurve>

- <https://hypertidy.github.io/bigcurve/>

- Report bugs at <https://github.com/hypertidy/bigcurve/issues>

## Author

**Maintainer**: Michael D. Sumner <mdsumner@gmail.com> \[copyright
holder\]

Authors:

- Michael D. Sumner <mdsumner@gmail.com> \[copyright holder\]

Other contributors:

- Simon Wotherspoon (wrote SGAT::trackMdpts, the ancestor of the
  spherical midpoint used here) \[contributor\]

- Mike Bostock (author of the 'd3-geo' adaptive resampling method ported
  here; see inst/COPYRIGHTS) \[contributor, copyright holder\]
