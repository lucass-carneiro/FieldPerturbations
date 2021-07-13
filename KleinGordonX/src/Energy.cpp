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
 * Energy.cpp
 * Compute the energy density of the wave equation.
 */

#include "KleinGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Energy;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);
  const GF3D2<CCTK_REAL> gf_epsilon(layout, epsilon);

  auto error_lambda =
      [&](const PointDesc &p) {
        CCTK_REAL dt_phi = gf_K_Phi(p.I);
        CCTK_REAL dx_phi =
            (gf_Phi(p.I + p.DI[0]) - gf_Phi(p.I - p.DI[0])) / (2 * p.dx);
        CCTK_REAL dy_phi =
            (gf_Phi(p.I + p.DI[1]) - gf_Phi(p.I - p.DI[1])) / (2 * p.dy);
        CCTK_REAL dz_phi =
            (gf_Phi(p.I + p.DI[2]) - gf_Phi(p.I - p.DI[2])) / (2 * p.dz);
        gf_epsilon(p.I) = (pow(dt_phi, 2) + pow(dx_phi, 2) + pow(dy_phi, 2) +
                       pow(dz_phi, 2)) /
                      2;
  };

  loop_int<1, 1, 1>(cctkGH, error_lambda);
}
