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
 * Error.cpp
 * Compute the error of the solution with respect to the analytic solution.
 */

#include "KleinGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  
  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);
  const GF3D2<CCTK_REAL> gf_Phi_err(layout, Phi_err);
  const GF3D2<CCTK_REAL> gf_K_Phi_err(layout, K_Phi_err);

  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    auto gaussian_lambda =
        [&](const PointDesc &p) {
          gf_Phi_err(p.I) = gf_Phi(p.I) - gaussian(t, p.x, p.y, p.z);
          gf_K_Phi_err(p.I) =
              gf_K_Phi(p.I) - timederiv(gaussian, dt)(t, p.x, p.y, p.z);
    };

    loop_all<1, 1, 1>(cctkGH, gaussian_lambda);
  }
}
