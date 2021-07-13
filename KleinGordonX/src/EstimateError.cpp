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
 * EstimateError.cpp
 * Computes the regrid error.
 */

#include "KlainGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<CCTK_REAL> gf_K_Phi(layout, K_Phi);
  const GF3D2<CCTK_REAL> gf_regrid_error(layout, regrid_error);

  auto regriderror_lambda =
      [&](const PointDesc &p) {
        CCTK_REAL base_phi = fabs(gf_Phi(p.I)) + fabs(phi_abs);
        CCTK_REAL errx_phi = fabs(gf_Phi(p.I - p.DI[0]) - 2 * gf_Phi(p.I) +
                                  gf_Phi(p.I + p.DI[0])) /
                             base_phi;
        CCTK_REAL erry_phi = fabs(gf_Phi(p.I - p.DI[1]) - 2 * gf_Phi(p.I) +
                                  gf_Phi(p.I + p.DI[1])) /
                             base_phi;
        CCTK_REAL errz_phi = fabs(gf_Phi(p.I - p.DI[2]) - 2 * gf_Phi(p.I) +
                                  gf_Phi(p.I + p.DI[2])) /
                             base_phi;
        CCTK_REAL base_psi = fabs(gf_K_Phi(p.I)) + fabs(psi_abs);
        CCTK_REAL errx_psi = fabs(gf_K_Phi(p.I - p.DI[0]) - 2 * gf_K_Phi(p.I) +
                                  gf_K_Phi(p.I + p.DI[0])) /
                             base_psi;
        CCTK_REAL erry_psi = fabs(gf_K_Phi(p.I - p.DI[1]) - 2 * gf_K_Phi(p.I) +
                                  gf_K_Phi(p.I + p.DI[1])) /
                             base_psi;
        CCTK_REAL errz_psi = fabs(gf_K_Phi(p.I - p.DI[2]) - 2 * gf_K_Phi(p.I) +
                                  gf_K_Phi(p.I + p.DI[2])) /
                             base_psi;
        gf_regrid_error(p.I) =
            errx_phi + erry_phi + errz_phi + errx_psi + erry_psi + errz_psi;
      }

  loop_int<1, 1, 1>(cctkGH, regriderror_lambda);
}
