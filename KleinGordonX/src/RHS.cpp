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
 * RHS.cpp
 * Compute the RHS of the wave equation
 */

#include "KleinGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_RHS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);

  const GF3D2<CCTK_REAL> gf_Phi_rhs(layout, Phi_rhs);
  const GF3D2<CCTK_REAL> gf_K_Phi_rhs(layout, K_Phi_rhs);

  auto rhs_lambda = [&](const PointDesc &p) {
    CCTK_REAL ddx_phi =
        (gf_Phi(p.I - p.DI[0]) - 2 * gf_Phi(p.I) + gf_Phi(p.I + p.DI[0])) /
        (p.dx * p.dx);
    CCTK_REAL ddy_phi =
        (gf_Phi(p.I - p.DI[1]) - 2 * gf_Phi(p.I) + gf_Phi(p.I + p.DI[1])) /
        (p.dx * p.dx);
    CCTK_REAL ddz_phi =
        (gf_Phi(p.I - p.DI[2]) - 2 * gf_Phi(p.I) + gf_Phi(p.I + p.DI[2])) /
        (p.dx * p.dx);
    gf_K_Phi_rhs(p.I) = ddx_phi + ddy_phi + ddz_phi -
                        pow(field_mass, 2) * gf_Phi(p.I) +
                        4 * M_PI * central_potential(t, p.x, p.y, p.z);
    gf_Phi_rhs(p.I) = gf_K_Phi(p.I);
  };

  loop_int<1, 1, 1>(cctkGH, rhs_lambda);
}
