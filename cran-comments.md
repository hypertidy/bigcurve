## Submission

This is a new release.

bigcurve densifies lines adaptively for a target map projection, following
the adaptive resampling method of 'd3-geo'. The part of the C++ code that
ports d3-geo's refinement test is credited in Authors@R (Mike Bostock, ctb
and cph) and its ISC licence notice is reproduced in inst/COPYRIGHTS.

The package links to the system 'PROJ' library (>= 8.0), found with
pkg-config by a configure script on Unix and macOS and through Rtools on
Windows.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Method references

The method is described by its author at <https://bost.ocks.org/mike/example/>
and implemented in 'd3-geo' <https://github.com/d3/d3-geo>; there is no
journal publication.
