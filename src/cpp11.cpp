// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// proj.cpp
void cpp_libproj_init_api();
extern "C" SEXP _bigcurve_cpp_libproj_init_api() {
  BEGIN_CPP11
    cpp_libproj_init_api();
    return R_NilValue;
  END_CPP11
}
// proj.cpp
list proj_coords(list xy, std::string from, std::string to);
extern "C" SEXP _bigcurve_proj_coords(SEXP xy, SEXP from, SEXP to) {
  BEGIN_CPP11
    return cpp11::as_sexp(proj_coords(cpp11::as_cpp<cpp11::decay_t<list>>(xy), cpp11::as_cpp<cpp11::decay_t<std::string>>(from), cpp11::as_cpp<cpp11::decay_t<std::string>>(to)));
  END_CPP11
}
// proj.cpp
doubles mid_pt1(std::vector<double> x1, std::vector<double> x2);
extern "C" SEXP _bigcurve_mid_pt1(SEXP x1, SEXP x2) {
  BEGIN_CPP11
    return cpp11::as_sexp(mid_pt1(cpp11::as_cpp<cpp11::decay_t<std::vector<double>>>(x1), cpp11::as_cpp<cpp11::decay_t<std::vector<double>>>(x2)));
  END_CPP11
}
// proj.cpp
list mid_pt_pairs_gc(doubles lon, doubles lat);
extern "C" SEXP _bigcurve_mid_pt_pairs_gc(SEXP lon, SEXP lat) {
  BEGIN_CPP11
    return cpp11::as_sexp(mid_pt_pairs_gc(cpp11::as_cpp<cpp11::decay_t<doubles>>(lon), cpp11::as_cpp<cpp11::decay_t<doubles>>(lat)));
  END_CPP11
}
// proj.cpp
list mid_pt_pairs_plane(doubles x, doubles y);
extern "C" SEXP _bigcurve_mid_pt_pairs_plane(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(mid_pt_pairs_plane(cpp11::as_cpp<cpp11::decay_t<doubles>>(x), cpp11::as_cpp<cpp11::decay_t<doubles>>(y)));
  END_CPP11
}
// proj.cpp
double gc_dist1(double lon1, double lat1, double lon2, double lat2);
extern "C" SEXP _bigcurve_gc_dist1(SEXP lon1, SEXP lat1, SEXP lon2, SEXP lat2) {
  BEGIN_CPP11
    return cpp11::as_sexp(gc_dist1(cpp11::as_cpp<cpp11::decay_t<double>>(lon1), cpp11::as_cpp<cpp11::decay_t<double>>(lat1), cpp11::as_cpp<cpp11::decay_t<double>>(lon2), cpp11::as_cpp<cpp11::decay_t<double>>(lat2)));
  END_CPP11
}
// proj.cpp
doubles dist_2_gc(doubles x0, doubles y0, doubles x1, doubles y1);
extern "C" SEXP _bigcurve_dist_2_gc(SEXP x0, SEXP y0, SEXP x1, SEXP y1) {
  BEGIN_CPP11
    return cpp11::as_sexp(dist_2_gc(cpp11::as_cpp<cpp11::decay_t<doubles>>(x0), cpp11::as_cpp<cpp11::decay_t<doubles>>(y0), cpp11::as_cpp<cpp11::decay_t<doubles>>(x1), cpp11::as_cpp<cpp11::decay_t<doubles>>(y1)));
  END_CPP11
}
// proj.cpp
list bisect_cpp(list xy, std::string from, std::string to);
extern "C" SEXP _bigcurve_bisect_cpp(SEXP xy, SEXP from, SEXP to) {
  BEGIN_CPP11
    return cpp11::as_sexp(bisect_cpp(cpp11::as_cpp<cpp11::decay_t<list>>(xy), cpp11::as_cpp<cpp11::decay_t<std::string>>(from), cpp11::as_cpp<cpp11::decay_t<std::string>>(to)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_bigcurve_bisect_cpp",           (DL_FUNC) &_bigcurve_bisect_cpp,           3},
    {"_bigcurve_cpp_libproj_init_api", (DL_FUNC) &_bigcurve_cpp_libproj_init_api, 0},
    {"_bigcurve_dist_2_gc",            (DL_FUNC) &_bigcurve_dist_2_gc,            4},
    {"_bigcurve_gc_dist1",             (DL_FUNC) &_bigcurve_gc_dist1,             4},
    {"_bigcurve_mid_pt1",              (DL_FUNC) &_bigcurve_mid_pt1,              2},
    {"_bigcurve_mid_pt_pairs_gc",      (DL_FUNC) &_bigcurve_mid_pt_pairs_gc,      2},
    {"_bigcurve_mid_pt_pairs_plane",   (DL_FUNC) &_bigcurve_mid_pt_pairs_plane,   2},
    {"_bigcurve_proj_coords",          (DL_FUNC) &_bigcurve_proj_coords,          3},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_bigcurve(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
