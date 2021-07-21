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
 * Boundaries.cpp
 * Implements outer boundary conditions for the ADM-decomposed scalar wave
 * equation and RHS boundary conditions
 */

#include "KleinGordonX.hpp"

using namespace Loop;

extern "C" void KleinGordonX::KleinGordonX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<CCTK_REAL> gf_K_Phi(layout, K_Phi);

  if (CCTK_EQUALS(bc_type, "NewRad")) {
    // Apply NewRad
  } else if (CCTK_EQUALS(bc_type, "none")) {
    // Do nothing
  } else if (CCTK_EQUALS(bc_type, "reflecting")) {

    auto dirichilet_lambda = [&](const PointDesc &p) {
      gf_Phi(p.I) = CCTK_REAL(0);
      gf_K_Phi(p.I) = CCTK_REAL(0);
    };

    loop_bnd<1, 1, 1>(cctkGH, dirichilet_lambda);
  }
}

extern "C" void KleinGordonX::KleinGordonX_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_RHSBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_Phi_rhs(layout, Phi_rhs);
  const GF3D2<CCTK_REAL> gf_K_Phi_rhs(layout, K_Phi_rhs);

  auto dirichilet_lambda = [&](const PointDesc &p) {
    gf_Phi_rhs(p.I) = CCTK_REAL(0);
    gf_K_Phi_rhs(p.I) = CCTK_REAL(0);
  };

  loop_bnd<1, 1, 1>(cctkGH, dirichilet_lambda);
}
