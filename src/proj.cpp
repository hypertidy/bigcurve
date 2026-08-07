#include <string>

#include <cpp11.hpp>
using namespace cpp11;

#include <proj.h>

// batch transform, used by rproj_xy() and the default tolerance heuristic
[[cpp11::register]]
list proj_coords(list xy, std::string from, std::string to) {
  doubles x = xy[0];
  doubles y = xy[1];

  PJ_CONTEXT* context = PJ_DEFAULT_CTX;

  PJ* raw = proj_create_crs_to_crs(context, from.c_str(), to.c_str(), NULL);
  if (raw == NULL) {
    int error_code = proj_context_errno(context);
    stop("Error creating transform: %s",
         proj_context_errno_string(context, error_code));
  }
  PJ* trans = proj_normalize_for_visualization(context, raw);
  if (trans != NULL) {
    proj_destroy(raw);
  } else {
    trans = raw;
  }

  writable::doubles xout(x);
  writable::doubles yout(y);
  size_t stride = sizeof(double);

  proj_trans_generic(
    trans, PJ_FWD,
    REAL(xout), stride, xout.size(),
    REAL(yout), stride, yout.size(),
    nullptr, stride, 0,
    nullptr, stride, 0
  );

  proj_destroy(trans);

  // do not error on per-point failures: out-of-domain points are normal
  // for this package; they come back non-finite and callers handle them

  writable::list out = {xout, yout};
  out.names() = {"x", "y"};
  return out;
}

// PROJ version string, handy for Sys/debug reporting
[[cpp11::register]]
std::string proj_version_cpp() {
  PJ_INFO info = proj_info();
  return std::string(info.release);
}
