#ifndef BIGCURVE_H
#define BIGCURVE_H

// bigcurve core: adaptive densification of segments under a nonlinear
// forward transform, after Bostock's d3-geo resampling.
//
// This header has no dependency on R or PROJ. The transform is injected
// via the Transform interface; bindings supply an adapter (see densify.cpp
// for the PROJ one). Keep it that way: this file is the piece that gets
// reused in Python (or rewritten in Rust) later.

#include <cmath>
#include <vector>

namespace bigcurve {

const double DEG2RAD = 0.017453292519943295;
const double RAD2DEG = 57.29577951308232;

// Forward transform lon/lat (degrees) -> planar x/y (projected units).
// Return false if the point cannot be transformed (out of domain etc.);
// implementations should also return false for non-finite output.
struct Transform {
  virtual bool fwd(double lon, double lat, double* x, double* y) = 0;
  virtual ~Transform() {}
};

// Spherical midpoint of the great circle arc between two lon/lat points,
// degrees in and out, longitude wrapped to [-180, 180).
inline void gc_midpoint(double lon1, double lat1, double lon2, double lat2,
                        double* lonm, double* latm) {
  double rlon1 = lon1 * DEG2RAD;
  double rlat1 = lat1 * DEG2RAD;
  double rlat2 = lat2 * DEG2RAD;
  double dlon = (lon2 - lon1) * DEG2RAD;
  double bx = std::cos(rlat2) * std::cos(dlon);
  double by = std::cos(rlat2) * std::sin(dlon);
  double cl1bx = std::cos(rlat1) + bx;
  *latm = std::atan2(std::sin(rlat1) + std::sin(rlat2),
                     std::sqrt(cl1bx * cl1bx + by * by)) * RAD2DEG;
  double lon = (rlon1 + std::atan2(by, cl1bx)) * RAD2DEG;
  // wrap that handles any input sign (fmod alone does not)
  *lonm = lon - 360.0 * std::floor((lon + 180.0) / 360.0);
}

// A vertex carries both its geographic and projected coordinates so that
// no point is ever transformed twice.
struct Vertex {
  double lon, lat;  // geographic, degrees
  double x, y;      // projected
  bool ok;          // projectable and finite
};

// cosine of the angular distance between two lon/lat points (degrees)
inline double cos_angular(double lon1, double lat1, double lon2, double lat2) {
  double rlat1 = lat1 * DEG2RAD;
  double rlat2 = lat2 * DEG2RAD;
  return std::sin(rlat1) * std::sin(rlat2) +
         std::cos(rlat1) * std::cos(rlat2) *
         std::cos((lon2 - lon1) * DEG2RAD);
}

class Resampler {
public:
  // tol: tolerance in projected (target) units. A segment is split when the
  // projected great-circle midpoint sits further than tol from the planar
  // chord, or lands parametrically far from the chord's middle (shear).
  Resampler(Transform& t, double tol, int max_depth)
    : t_(t), tol2_(tol * tol), max_depth_(max_depth) {}

  Vertex project(double lon, double lat) {
    Vertex v;
    v.lon = lon;
    v.lat = lat;
    v.ok = t_.fwd(lon, lat, &v.x, &v.y) &&
           std::isfinite(v.x) && std::isfinite(v.y);
    return v;
  }

  // Densify one segment a-b, appending interior vertices (in order, not
  // including the endpoints) to out. Segments with an unprojectable
  // endpoint are left alone; unprojectable midpoints stop refinement of
  // that branch rather than erroring.
  void edge(const Vertex& a, const Vertex& b, std::vector<Vertex>& out) {
    if (!a.ok || !b.ok) return;
    split(a, b, max_depth_, out);
  }

private:
  void split(const Vertex& a, const Vertex& b, int depth,
             std::vector<Vertex>& out) {
    if (depth <= 0) return;

    // arcs wider than ~30 degrees always split: the midpoint metric below
    // is blind to symmetric curves whose midpoint lands on the chord
    // (e.g. an oblique great circle through its own inflection), so it is
    // only allowed to adjudicate once arcs are short (d3's cosMinDistance)
    const double COS_MIN_DISTANCE = 0.8660254037844387;  // cos(30 deg)
    bool wide = cos_angular(a.lon, a.lat, b.lon, b.lat) < COS_MIN_DISTANCE;

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double d2 = dx * dx + dy * dy;
    // chord shorter than 2 * tol and arc not wide: nothing left to resolve
    if (!wide && d2 <= 4.0 * tol2_) return;

    double lonm, latm;
    gc_midpoint(a.lon, a.lat, b.lon, b.lat, &lonm, &latm);
    Vertex m = project(lonm, latm);
    if (!m.ok) return;

    bool go = wide;
    if (!go && d2 > 0.0) {
      double dx2 = m.x - a.x;
      double dy2 = m.y - a.y;
      double dz = dx * dy2 - dy * dx2;           // 2 * signed area
      double perp2 = (dz * dz) / d2;             // squared perpendicular dist
      double tpar = (dx * dx2 + dy * dy2) / d2;  // parametric position
      go = perp2 > tol2_ || tpar < 0.2 || tpar > 0.8;
    }

    if (go) {
      split(a, m, depth - 1, out);
      out.push_back(m);
      split(m, b, depth - 1, out);
    }
  }

  Transform& t_;
  double tol2_;
  int max_depth_;
};

}  // namespace bigcurve

#endif
