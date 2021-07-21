#include "KleinGordonX.hpp"

using namespace Loop;

CCTK_REAL KleinGordonX::central_potential(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                          CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  return -1.0 / pow(1.0, dim) * spline_potential(r / 1.0);
}
