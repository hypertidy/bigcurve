#include <string>
#include <vector>

#include <cpp11.hpp>

#include <proj.h>

#include "bigcurve.h"

using namespace cpp11;

namespace {

// RAII PROJ adapter: one normalized PJ for the lifetime of a densify call.
class ProjTransform : public bigcurve::Transform {
public:
  ProjTransform(const std::string& from, const std::string& to) {
    PJ_CONTEXT* ctx = PJ_DEFAULT_CTX;
    PJ* raw = proj_create_crs_to_crs(ctx, from.c_str(), to.c_str(), NULL);
    if (raw == NULL) {
      stop("Error creating transform: %s",
           proj_context_errno_string(ctx, proj_context_errno(ctx)));
    }
    // always lon/lat in, easting/northing out, regardless of authority
    // axis order (so EPSG:4326 behaves like OGC:CRS84 here)
    PJ* norm = proj_normalize_for_visualization(ctx, raw);
    if (norm != NULL) {
      proj_destroy(raw);
      trans_ = norm;
    } else {
      trans_ = raw;
    }
  }

  ~ProjTransform() { proj_destroy(trans_); }

  bool fwd(double lon, double lat, double* x, double* y) override {
    PJ_COORD in;
    in.xyzt.x = lon;
    in.xyzt.y = lat;
    in.xyzt.z = 0.0;
    in.xyzt.t = HUGE_VAL;
    PJ_COORD out = proj_trans(trans_, PJ_FWD, in);
    *x = out.xy.x;
    *y = out.xy.y;
    return std::isfinite(*x) && std::isfinite(*y);
  }

private:
  PJ* trans_;
};

writable::doubles to_doubles(const std::vector<double>& v) {
  writable::doubles out((R_xlen_t)v.size());
  for (R_xlen_t i = 0; i < (R_xlen_t)v.size(); i++) out[i] = v[(size_t)i];
  return out;
}

writable::integers to_integers(const std::vector<int>& v) {
  writable::integers out((R_xlen_t)v.size());
  for (R_xlen_t i = 0; i < (R_xlen_t)v.size(); i++) out[i] = v[(size_t)i];
  return out;
}

}  // namespace

// Densify an ordered path of lon/lat vertices. Returns list(x, y) of the
// densified path in geographic coordinates, endpoints preserved.
[[cpp11::register]]
list densify_path_cpp(doubles lon, doubles lat,
                      std::string from, std::string to,
                      double tol, int max_depth) {
  ProjTransform t(from, to);
  bigcurve::Resampler rs(t, tol, max_depth);

  R_xlen_t n = lon.size();
  std::vector<bigcurve::Vertex> v((size_t)n);
  for (R_xlen_t i = 0; i < n; i++) {
    v[(size_t)i] = rs.project(lon[i], lat[i]);
  }

  std::vector<double> ox, oy;
  ox.reserve((size_t)n * 2);
  oy.reserve((size_t)n * 2);
  std::vector<bigcurve::Vertex> interior;

  for (R_xlen_t i = 0; i + 1 < n; i++) {
    ox.push_back(v[(size_t)i].lon);
    oy.push_back(v[(size_t)i].lat);
    interior.clear();
    rs.edge(v[(size_t)i], v[(size_t)i + 1], interior);
    for (size_t k = 0; k < interior.size(); k++) {
      ox.push_back(interior[k].lon);
      oy.push_back(interior[k].lat);
    }
  }
  if (n > 0) {
    ox.push_back(v[(size_t)(n - 1)].lon);
    oy.push_back(v[(size_t)(n - 1)].lat);
  }

  writable::list out({to_doubles(ox), to_doubles(oy)});
  out.names() = {"x", "y"};
  return out;
}

// Densify a segment mesh: vertices (x, y) and 0-based edge index pairs
// (s0, s1). Returns list(x, y, s0, s1, parent) with new vertices appended,
// each input edge replaced by its refined chain (0-based indices), and
// parent giving the 0-based input edge each output edge descends from --
// the hook for carrying per-segment attributes (feature, ring, arc)
// through refinement.
[[cpp11::register]]
list densify_mesh_cpp(doubles x, doubles y,
                      integers s0, integers s1,
                      std::string from, std::string to,
                      double tol, int max_depth) {
  ProjTransform t(from, to);
  bigcurve::Resampler rs(t, tol, max_depth);

  R_xlen_t nv = x.size();
  std::vector<bigcurve::Vertex> v((size_t)nv);
  for (R_xlen_t i = 0; i < nv; i++) {
    v[(size_t)i] = rs.project(x[i], y[i]);
  }

  std::vector<double> vx(x.begin(), x.end());
  std::vector<double> vy(y.begin(), y.end());
  std::vector<int> e0, e1, parent;
  e0.reserve((size_t)s0.size() * 2);
  e1.reserve((size_t)s0.size() * 2);
  parent.reserve((size_t)s0.size() * 2);
  std::vector<bigcurve::Vertex> interior;

  for (R_xlen_t k = 0; k < s0.size(); k++) {
    int a = s0[k];
    int b = s1[k];
    if (a < 0 || b < 0 || a >= (int)nv || b >= (int)nv) {
      stop("edge index out of range");
    }
    interior.clear();
    rs.edge(v[(size_t)a], v[(size_t)b], interior);
    int prev = a;
    for (size_t j = 0; j < interior.size(); j++) {
      vx.push_back(interior[j].lon);
      vy.push_back(interior[j].lat);
      int idx = (int)vx.size() - 1;
      e0.push_back(prev);
      e1.push_back(idx);
      parent.push_back((int)k);
      prev = idx;
    }
    e0.push_back(prev);
    e1.push_back(b);
    parent.push_back((int)k);
  }

  writable::list out({to_doubles(vx), to_doubles(vy),
                      to_integers(e0), to_integers(e1),
                      to_integers(parent)});
  out.names() = {"x", "y", "s0", "s1", "parent"};
  return out;
}
