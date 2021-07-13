/*
 *  KleinGordonX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordonX.
 *
 *  KleinGordonX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordonX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KleinGordonX.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Initialize.cpp
 *  Fill the grid functions with initial data.
 */

#include "KleinGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<CCTK_REAL> gf_K_Phi(layout, K_Phi);

  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    auto gaussian_lambda =
        [&](const PointDesc &p) {
          gf_Phi(p.I) = gaussian(p.x, p.y, p.z);
          gf_K_Phi(p.I) = CCTK_REAL(0);
        }

    loop_int<1, 1, 1>(cctkGH, gaussian_lambda);
  }
}

CCTK_REAL KleinGordonX::gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                 CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  // u(t,r) = (f(r-t) - f(r+t)) / r
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  auto f = [&](auto x) { return exp(-0.5 * pow(x / width, 2)); };
  auto fx = [&](auto x) { return -x / pow(width, 2) * f(x); };

  // Use L'Hôpital's rule for small r
  if (r < 1.0e-8)
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}
