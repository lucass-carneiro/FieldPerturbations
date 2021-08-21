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
#include "Derivative.hpp"

namespace KleinGordonX {

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_int;
using Loop::PointDesc;

extern "C" void KleinGordonX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Energy;
  DECLARE_CCTK_PARAMETERS;

  const vect<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);
  const GF3D2<CCTK_REAL> gf_epsilon(layout, epsilon);

  const Derivative<Stencil::c3> d_Phi(gf_Phi);

  auto error_lambda = [&](const PointDesc &p) {
    auto grad_Phi = d_Phi.grad(p);

    const CCTK_REAL dt_Phi = gf_K_Phi(p.I);
    const CCTK_REAL dx_Phi = grad_Phi[0];
    const CCTK_REAL dy_Phi = grad_Phi[1];
    const CCTK_REAL dz_Phi = grad_Phi[2];

    gf_epsilon(p.I) =
        pow2(dt_Phi) + pow2(dx_Phi) + pow2(dy_Phi) + pow2(dz_Phi) / 2;
  };

  loop_int<0, 0, 0>(cctkGH, error_lambda);
}
} // namespace KleinGordonX
