// check.cpp -- the C++ engine against the canon tables of bigcurve-rs.
//
// hypertidy/bigcurve-rs is the canonical bigcurve engine; src/bigcurve.h
// (frozen, bug fixes only) follows it. canon/canon-tables.txt is written
// by bigcurve-rs (tests/canon.rs): inputs, settings, every forward call
// the Rust engine made, and its output. This program replays those calls
// through the header, driven exactly as src/densify.cpp drives it, and
// requires the same output: bitwise on Linux with glibc (where the tables
// are written), elsewhere the same vertex counts and segments with
// coordinates to 1e-9 (the engine's trigonometry comes from the
// platform's maths library). Cases with Planar edges or ellipsoid
// placement are Rust-only for now and are counted, not checked.
//
// Build and run from the package root (see canon/README.md):
//
//     c++ -std=c++17 -O2 -o canon/check canon/check.cpp && ./canon/check

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../src/bigcurve.h"

namespace {

#if defined(__linux__) && defined(__GLIBC__)
const bool STRICT = true;
#else
const bool STRICT = false;
#endif
const double CLOSE = 1e-9;

double unhex(const std::string& s) {
  std::uint64_t u = std::stoull(s, nullptr, 16);
  double v;
  std::memcpy(&v, &u, sizeof v);
  return v;
}

std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return u;
}

struct Call {
  bool ok;
  double x, y;
};

struct Case {
  std::string name, kind, edge, place;
  double tol = 0;
  int max_depth = 16;
  std::vector<double> x, y;
  std::vector<int> s0, s1;
  std::map<std::pair<std::uint64_t, std::uint64_t>, Call> calls;
  std::vector<double> ox, oy;
  std::vector<int> o0, o1, op;
};

class Replay : public bigcurve::Transform {
public:
  Replay(const Case& c) : c_(c) {}
  bool fwd(double lon, double lat, double* x, double* y) override {
    auto it = c_.calls.find({bits(lon), bits(lat)});
    if (it == c_.calls.end() && !STRICT) {
      for (auto j = c_.calls.begin(); j != c_.calls.end(); ++j) {
        double a, b;
        std::uint64_t u = j->first.first, v = j->first.second;
        std::memcpy(&a, &u, sizeof a);
        std::memcpy(&b, &v, sizeof b);
        if (std::fabs(a - lon) < CLOSE && std::fabs(b - lat) < CLOSE) {
          it = j;
          break;
        }
      }
    }
    if (it == c_.calls.end()) {
      missing = true;
      *x = NAN;
      *y = NAN;
      return false;
    }
    *x = it->second.x;
    *y = it->second.y;
    return it->second.ok;
  }
  bool missing = false;

private:
  const Case& c_;
};

struct Out {
  std::vector<double> x, y;
  std::vector<int> s0, s1, parent;
};

// densify_path_cpp in src/densify.cpp, without R
Out path(bigcurve::Transform& t, const Case& c) {
  bigcurve::Resampler rs(t, c.tol, c.max_depth);
  std::size_t n = c.x.size();
  std::vector<bigcurve::Vertex> v(n);
  for (std::size_t i = 0; i < n; i++) v[i] = rs.project(c.x[i], c.y[i]);
  Out o;
  std::vector<bigcurve::Vertex> interior;
  for (std::size_t i = 0; i + 1 < n; i++) {
    o.x.push_back(v[i].lon);
    o.y.push_back(v[i].lat);
    interior.clear();
    rs.edge(v[i], v[i + 1], interior);
    for (auto& w : interior) {
      o.x.push_back(w.lon);
      o.y.push_back(w.lat);
    }
  }
  if (n > 0) {
    o.x.push_back(v[n - 1].lon);
    o.y.push_back(v[n - 1].lat);
  }
  return o;
}

// densify_mesh_cpp in src/densify.cpp, without R
Out mesh(bigcurve::Transform& t, const Case& c) {
  bigcurve::Resampler rs(t, c.tol, c.max_depth);
  int nv = (int)c.x.size();
  std::vector<bigcurve::Vertex> v((std::size_t)nv);
  for (int i = 0; i < nv; i++) v[(std::size_t)i] = rs.project(c.x[(std::size_t)i], c.y[(std::size_t)i]);
  Out o;
  o.x = c.x;
  o.y = c.y;
  std::vector<bigcurve::Vertex> interior;
  struct Seen { int from; int first; int n; };
  std::unordered_map<long long, Seen> seen;
  std::vector<int> chain;
  for (std::size_t k = 0; k < c.s0.size(); k++) {
    int a = c.s0[k];
    int b = c.s1[k];
    long long key = (long long)std::min(a, b) * (long long)nv + (long long)std::max(a, b);
    chain.clear();
    auto hit = seen.find(key);
    if (hit != seen.end()) {
      const Seen& h = hit->second;
      for (int j = 0; j < h.n; j++) chain.push_back(h.first + j);
      if (h.from != a) std::reverse(chain.begin(), chain.end());
    } else {
      interior.clear();
      rs.edge(v[(std::size_t)a], v[(std::size_t)b], interior);
      int first = (int)o.x.size();
      for (std::size_t j = 0; j < interior.size(); j++) {
        o.x.push_back(interior[j].lon);
        o.y.push_back(interior[j].lat);
        chain.push_back(first + (int)j);
      }
      seen.emplace(key, Seen{a, first, (int)interior.size()});
    }
    int prev = a;
    for (int idx : chain) {
      o.s0.push_back(prev);
      o.s1.push_back(idx);
      o.parent.push_back((int)k);
      prev = idx;
    }
    o.s0.push_back(prev);
    o.s1.push_back(b);
    o.parent.push_back((int)k);
  }
  return o;
}

std::vector<Case> read(const std::string& file) {
  std::ifstream in(file);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", file.c_str());
    std::exit(2);
  }
  std::vector<Case> out;
  Case c;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream w(line);
    std::string tag;
    if (!(w >> tag) || tag[0] == '#') continue;
    if (tag == "case") {
      c = Case();
      w >> c.name;
    } else if (tag == "kind") {
      w >> c.kind;
    } else if (tag == "edge") {
      w >> c.edge;
    } else if (tag == "place") {
      w >> c.place;
    } else if (tag == "tolerance") {
      std::string h;
      w >> h;
      c.tol = unhex(h);
    } else if (tag == "max_depth") {
      w >> c.max_depth;
    } else if (tag == "in") {
      std::string a, b;
      w >> a >> b;
      c.x.push_back(unhex(a));
      c.y.push_back(unhex(b));
    } else if (tag == "seg") {
      int a, b;
      w >> a >> b;
      c.s0.push_back(a);
      c.s1.push_back(b);
    } else if (tag == "fwd") {
      std::string a, b, ok, px, py;
      w >> a >> b >> ok >> px >> py;
      c.calls[{std::stoull(a, nullptr, 16), std::stoull(b, nullptr, 16)}] =
          Call{ok == "ok", unhex(px), unhex(py)};
    } else if (tag == "out") {
      std::string a, b;
      w >> a >> b;
      c.ox.push_back(unhex(a));
      c.oy.push_back(unhex(b));
    } else if (tag == "oseg") {
      int a, b, p;
      w >> a >> b >> p;
      c.o0.push_back(a);
      c.o1.push_back(b);
      c.op.push_back(p);
    } else if (tag == "end") {
      out.push_back(c);
    } else {
      std::fprintf(stderr, "unknown record %s\n", tag.c_str());
      std::exit(2);
    }
  }
  return out;
}

bool same(double a, double b) {
  return STRICT ? bits(a) == bits(b) : std::fabs(a - b) < CLOSE;
}

// returns an empty string when the case agrees
std::string check(const Case& c) {
  Replay t(c);
  Out o = c.kind == "mesh" ? mesh(t, c) : path(t, c);
  if (t.missing) return "asked for a point the Rust engine never moved";
  if (o.x.size() != c.ox.size()) {
    return "vertex count " + std::to_string(o.x.size()) + " != " + std::to_string(c.ox.size());
  }
  for (std::size_t i = 0; i < o.x.size(); i++) {
    if (!same(o.x[i], c.ox[i]) || !same(o.y[i], c.oy[i])) {
      char buf[160];
      std::snprintf(buf, sizeof buf, "vertex %zu: (%.17g, %.17g) != (%.17g, %.17g)", i, o.x[i],
                    o.y[i], c.ox[i], c.oy[i]);
      return buf;
    }
  }
  if (o.s0 != c.o0 || o.s1 != c.o1 || o.parent != c.op) return "segments or parents differ";
  return "";
}

}  // namespace

int main(int argc, char** argv) {
  std::string file = argc > 1 ? argv[1] : "canon/canon-tables.txt";
  std::vector<Case> cases = read(file);
  int checked = 0, skipped = 0, failed = 0;
  for (const Case& c : cases) {
    if (c.edge != "greatcircle" || c.place != "sphere") {
      skipped++;
      continue;
    }
    checked++;
    std::string why = check(c);
    if (!why.empty()) {
      failed++;
      std::printf("FAIL %s: %s\n", c.name.c_str(), why.c_str());
    }
  }
  std::printf("canon tables: %d checked, %d failed, %d Rust-only (%s)\n", checked, failed, skipped,
              STRICT ? "bitwise" : "to 1e-9");
  if (checked == 0) return 1;
  return failed == 0 ? 0 : 1;
}
