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

#include "KleinGordonX.hpp"

namespace KleinGordonX {

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_int;
using Loop::PointDesc;

extern "C" void KleinGordonX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  /* Note that CarpetX's regrid_error grid function is cell centered.
   * Because our own variables are vertex centered, we must define two
   * separate grid layout objects: One cell centered and one vertex centered
   * and pass each to the corresponding GF.
   */
  const vect<int, dim> indextypecc = {1, 1, 1};
  const vect<int, dim> indextype = {0, 0, 0};

  const GF3D2layout layoutcc(cctkGH, indextypecc);
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);
  const GF3D2<CCTK_REAL> gf_regrid_error(layoutcc, regrid_error);

  auto regriderror_lambda = [&](const PointDesc &p) {
    const CCTK_REAL base_Phi = fabs(gf_Phi(p.I)) + fabs(Phi_abs);
    const CCTK_REAL errx_Phi =
        fabs(gf_Phi(p.I + p.DI[0]) - gf_Phi(p.I)) / base_Phi;
    const CCTK_REAL erry_Phi =
        fabs(gf_Phi(p.I + p.DI[1]) - gf_Phi(p.I)) / base_Phi;
    const CCTK_REAL errz_Phi =
        fabs(gf_Phi(p.I + p.DI[2]) - gf_Phi(p.I)) / base_Phi;

    const CCTK_REAL base_K_Phi = fabs(gf_K_Phi(p.I)) + fabs(K_Phi_abs);
    const CCTK_REAL errx_K_Phi =
        fabs(gf_K_Phi(p.I + p.DI[0]) - gf_K_Phi(p.I)) / base_K_Phi;
    const CCTK_REAL erry_K_Phi =
        fabs(gf_K_Phi(p.I + p.DI[1]) - gf_K_Phi(p.I)) / base_K_Phi;
    const CCTK_REAL errz_K_Phi =
        fabs(gf_K_Phi(p.I + p.DI[2]) - gf_K_Phi(p.I)) / base_K_Phi;

    gf_regrid_error(p.I) =
        errx_Phi + erry_Phi + errz_Phi + errx_K_Phi + erry_K_Phi + errz_K_Phi;
  };

  loop_int<1, 1, 1>(cctkGH, regriderror_lambda);
}

}
