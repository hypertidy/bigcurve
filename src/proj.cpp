#include <vector>
#include <math.h>


#include <cpp11.hpp>
using namespace cpp11;

// needed in every file that uses proj_*() functions
#include "libproj.h"

// needed exactly once in your package or Rcpp script
// contains all the function pointers and the
// implementation of the function to initialize them
// (`libproj_init_api()`)
#include "libproj.c"

// this function needs to be called once before any proj_*() functions
// are called (e.g., in .onLoad() for your package)
[[cpp11::linking_to(libproj)]]
[[cpp11::register]]
void cpp_libproj_init_api() {
  libproj_init_api();
}

// regular C or C++ code that uses proj_()* functions!
[[cpp11::register]]
list proj_coords(list xy, std::string from, std::string to) {
  doubles x = xy[0];
  doubles y = xy[1];

  PJ_CONTEXT* context = PJ_DEFAULT_CTX;

  PJ* trans = proj_create_crs_to_crs(context, from.c_str(), to.c_str(), NULL);
  if (trans == NULL) {
    int error_code = proj_context_errno(context);
    stop("Error creating transform: %s", proj_context_errno_string(context, error_code));
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

  int error_code = proj_errno(trans);
  proj_destroy(trans);

  if (error_code != 0) {
    stop("Error transforming coords: %s", proj_context_errno_string(context, error_code));
  }

  writable::list out = {xout, yout};
  out.names() = {"x", "y"};
  return out;
}


# define M_PI           3.14159265358979323846  /* pi */
# define M_RAD          6378137.0

// std::vector<double> geocentric_1(double lon, double lat) {
//
//   double cosLat = cos(lat * M_PI/180);
//   double sinLat = sin(lat * M_PI/180);
//   double cosLon = cos(lon * M_PI/180);
//   double sinLon = sin(lon * M_PI/180);
//
//   double x = M_RAD * cosLat * cosLon;
//   double y = M_RAD * cosLat * sinLon;
//   double z = M_RAD * sinLat;
//
//  std::vector<double> vout = {x, y, z};
//  return vout;
// }
// [[cpp11::register]]
// list geocentric_cpp(list xy, double rad, double exag)
// {
//   doubles x = xy[0];
//   doubles y = xy[1];
//
//   writable::doubles xout(x.size());
//   writable::doubles yout(x.size());
//   writable::doubles zout(x.size());
//
//   std::vector<double> xyz(3);
//   for (int i = 0; i < x.size();i++) {
//     xyz = geocentric_1(x[i], y[i]);
//     xout[i] = xyz[0];
//     yout[i] = xyz[1];
//     zout[i] = xyz[2];
//   }
//   writable::list out = {xout, yout, zout};
//   out.names() = {"x", "y", "z"};
//   return out;
// }
//
//
// double clamp1(double x) {
//   if (x > 1.0) {
//     x = 1.0;
//   }
//   if (x < -1.0) {
//     x = -1.0;
//   }
//   return x;
// }
//
// std::vector<double> longlat1(std::vector<double> x) {
//   double Zr = clamp1(x[2]/M_RAD);
//   std::vector<double> xout(2);
//   xout[0] = atan2(x[1], x[0]) * 180.0/M_PI;
//   xout[1] = asin(Zr) * 180.0/M_PI;
//   return xout;
// }

[[cpp11::register]]
doubles mid_pt1(std::vector<double> x1, std::vector<double> x2) {
  double dlon = x2[0] * M_PI/180.0 - x1[0]* M_PI/180.0;
  double lon1 = x1[0]* M_PI/180.0;
  double lat1 = x1[1]* M_PI/180.0;
  double lat2 = x2[1]* M_PI/180.0;
  double bx = cos(lat2) * cos(dlon);
  double by =  cos(lat2) * sin(dlon);
  double lat = atan2(sin(lat1) + sin(lat2), sqrt(pow(cos(lat1) + bx, 2.0) + by * by)) * 180.0 / M_PI;
  double lon =  (lon1 + atan2(by, cos(lat1) + bx)) * 180.0 / M_PI;

  writable::doubles out(2);
  out[0] = lon;
  out[1] = lat;
  return out;
}
[[cpp11::register]]
double gc_dist1(double lon1, double lat1, double lon2, double lat2) {

  double rad =  M_PI/180.0;

  double s1 = sin(rad * lat1);
  double s2 = sin(rad * lat2);
  double pin = sqrt(1 - s1 * s1) * sqrt(1 - s2 * s2) * cos(rad * (lon2 - lon1)) + s1 * s2;
  if (pin > 1.0) {
    pin = 1.0;
  }
  return M_RAD * acos(pin);
}




// regular C or C++ code that uses proj_()* functions!
[[cpp11::register]]
doubles bisect1(list xy, std::string from, std::string to) {
  doubles x = xy[0];
  doubles y = xy[1];

  PJ_CONTEXT* context = PJ_DEFAULT_CTX;

  PJ* trans = proj_create_crs_to_crs(context, from.c_str(), to.c_str(), NULL);

  if (trans == NULL) {
    int error_code = proj_context_errno(context);
    stop("Error creating transform: %s", proj_context_errno_string(context, error_code));
  }

  writable::doubles xout(x);
  writable::doubles yout(y);
  size_t stride = sizeof(double);

  // get the arc centroid
  doubles arc_cent = mid_pt1({x[0], y[0]}, {x[1], y[1]});
 // get the projected centroid
  proj_trans_generic(
    trans, PJ_FWD,
    REAL(xout), stride, xout.size(),
    REAL(yout), stride, yout.size(),
    nullptr, stride, 0,
    nullptr, stride, 0
  );

  writable::doubles xcent(1);
  writable::doubles ycent(1);

  xcent[0] = (xout[0] + xout[1])/2.0;
  ycent[0] = (yout[0] + yout[1])/2.0;

  proj_trans_generic(
    trans, PJ_INV,
    REAL(xcent), stride, 1,
    REAL(ycent), stride, 1,
    nullptr, stride, 0,
    nullptr, stride, 0
  );


  double dist = gc_dist1(arc_cent[0], arc_cent[1], xcent[0], ycent[1]);


  int error_code = proj_errno(trans);
  proj_destroy(trans);

  if (error_code != 0) {
    stop("Error transforming coords: %s", proj_context_errno_string(context, error_code));
  }

  writable::doubles out(5);
  out[0] = arc_cent[0];
  out[1] = arc_cent[1];
  out[2] = (double)xcent[0];
  out[3] = (double)ycent[0];
  out[4] = dist;

  //out.names() = {"dist"};
  return out;
}
