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
 * CalcTmunu.cpp
 * Compute the energy momentum tensor of the wave equation.
 */

#include "Derivative.hpp"
#include "KleinGordonX.hpp"

#include <iostream>

namespace KleinGordonX {

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_int;
using Loop::PointDesc;

extern "C" void KleinGordonX_CalcTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_CalcTmunu;
  DECLARE_CCTK_PARAMETERS;

  // Grid layout -----------------------------------------------
  const vect<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  // Field grid functions --------------------------------------
  const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);

  // ADM variables ---------------------------------------------
  const GF3D2<const CCTK_REAL> gf_alp(layout, alp);

  const GF3D2<const CCTK_REAL> gf_betax(layout, betax);
  const GF3D2<const CCTK_REAL> gf_betay(layout, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(layout, betaz);

  const GF3D2<const CCTK_REAL> gf_gxx(layout, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout, gzz);

  // Tmunu variables -------------------------------------------
  const GF3D2<CCTK_REAL> gf_Ttt(layout, eTtt);
  const GF3D2<CCTK_REAL> gf_Ttx(layout, eTtx);
  const GF3D2<CCTK_REAL> gf_Tty(layout, eTty);
  const GF3D2<CCTK_REAL> gf_Ttz(layout, eTtz);
  const GF3D2<CCTK_REAL> gf_Txx(layout, eTxx);
  const GF3D2<CCTK_REAL> gf_Txy(layout, eTxy);
  const GF3D2<CCTK_REAL> gf_Txz(layout, eTxz);
  const GF3D2<CCTK_REAL> gf_Tyy(layout, eTyy);
  const GF3D2<CCTK_REAL> gf_Tyz(layout, eTyz);
  const GF3D2<CCTK_REAL> gf_Tzz(layout, eTzz);

  // Derivative objects -----------------------------------------
  const Derivative<Stencil::c3> d_Phi(gf_Phi);

  auto Tmunu_lambda = [&](const PointDesc &p) {
    // Local ADM variables
    const CCTK_REAL alpL = gf_alp(p.I);

    const CCTK_REAL betaxL = gf_betax(p.I);
    const CCTK_REAL betayL = gf_betay(p.I);
    const CCTK_REAL betazL = gf_betaz(p.I);

    const CCTK_REAL gxxL = gf_gxx(p.I);
    const CCTK_REAL gxyL = gf_gxy(p.I);
    const CCTK_REAL gxzL = gf_gxz(p.I);
    const CCTK_REAL gyyL = gf_gyy(p.I);
    const CCTK_REAL gyzL = gf_gyz(p.I);
    const CCTK_REAL gzzL = gf_gzz(p.I);

    // Local grid functions
    const CCTK_REAL PhiL = gf_Phi(p.I);
    const CCTK_REAL K_PhiL = gf_K_Phi(p.I);

    // Computing the inverse metric
    const CCTK_REAL gdetL = -(gxzL * gxzL * gyyL) + 2 * gxyL * gxzL * gyzL -
                            gxxL * gyzL * gyzL - gxyL * gxyL * gzzL +
                            gxxL * gyyL * gzzL;
    const CCTK_REAL igxxL = (-gyzL * gyzL + gyyL * gzzL) / gdetL;
    const CCTK_REAL igxyL = (gxzL * gyzL - gxyL * gzzL) / gdetL;
    const CCTK_REAL igxzL = (-(gxzL * gyyL) + gxyL * gyzL) / gdetL;
    const CCTK_REAL igyyL = (-gxzL * gxzL + gxxL * gzzL) / gdetL;
    const CCTK_REAL igyzL = (gxyL * gxzL - gxxL * gyzL) / gdetL;
    const CCTK_REAL igzzL = (-gxyL * gxyL + gxxL * gyyL) / gdetL;

    // Derivatives of Phi
    auto grad_Phi = d_Phi.grad(p);

    const CCTK_REAL d_x_Phi = grad_Phi[0];
    const CCTK_REAL d_y_Phi = grad_Phi[1];
    const CCTK_REAL d_z_Phi = grad_Phi[2];
    const CCTK_REAL d_t_Phi = (betaxL * d_x_Phi + betayL * d_y_Phi + betazL * d_z_Phi) - 2 * alpL * K_PhiL;

    // 4 metric components (covariant)
    const CCTK_REAL gtt_4 =
        -alpL * alpL + gxxL * betaxL * betaxL + gxyL * betaxL * betayL * 2.0 +
        gxzL * betaxL * betazL * 2.0 + gyyL * betayL * betayL +
        gyzL * betayL * betazL * 2.0 + gzzL * betazL * betazL;
    const CCTK_REAL gtx_4 = gxxL * betaxL + gxyL * betayL + gxzL * betazL;
    const CCTK_REAL gty_4 = gxyL * betaxL + gyyL * betayL + gyzL * betazL;
    const CCTK_REAL gtz_4 = gxzL * betaxL + gyzL * betayL + gzzL * betazL;
    const CCTK_REAL gxx_4 = gxxL;
    const CCTK_REAL gxy_4 = gxyL;
    const CCTK_REAL gxz_4 = gxzL;
    const CCTK_REAL gyy_4 = gyyL;
    const CCTK_REAL gyz_4 = gyzL;
    const CCTK_REAL gzz_4 = gzzL;

    // 4 metric components (contravariant)
    const CCTK_REAL igtt_4 = -1.0 / (alpL * alpL);
    const CCTK_REAL igtx_4 = -igtt_4 * betaxL;
    const CCTK_REAL igty_4 = -igtt_4 * betayL;
    const CCTK_REAL igtz_4 = -igtt_4 * betazL;
    const CCTK_REAL igxx_4 = igxxL + igtt_4 * betaxL * betaxL;
    const CCTK_REAL igxy_4 = igxyL + igtt_4 * betaxL * betayL;
    const CCTK_REAL igxz_4 = igxzL + igtt_4 * betaxL * betazL;
    const CCTK_REAL igyy_4 = igyyL + igtt_4 * betayL * betayL;
    const CCTK_REAL igyz_4 = igyzL + igtt_4 * betayL * betazL;
    const CCTK_REAL igzz_4 = igzzL + igtt_4 * betazL * betazL;

    // Inner product of derivatives
    const CCTK_REAL contracted_part =
        igtt_4 * d_t_Phi * d_t_Phi + igtx_4 * d_t_Phi * d_x_Phi * 2 +
        igty_4 * d_t_Phi * d_y_Phi * 2 + igtz_4 * d_t_Phi * d_z_Phi * 2 +
        igxx_4 * d_x_Phi * d_x_Phi + igxy_4 * d_x_Phi * d_y_Phi * 2 +
        igxz_4 * d_x_Phi * d_z_Phi * 2 + igyy_4 * d_y_Phi * d_y_Phi +
        igyz_4 * d_y_Phi * d_z_Phi * 2 + igzz_4 * d_z_Phi * d_z_Phi;

    // T_{\mu \nu} components
    gf_Ttt(p.I) += d_t_Phi * d_t_Phi - 0.5 * gtx_4 * contracted_part + 0.5 * gtx_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Ttx(p.I) += d_t_Phi * d_x_Phi - 0.5 * gtx_4 * contracted_part + 0.5 * gtx_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Tty(p.I) += d_t_Phi * d_y_Phi - 0.5 * gty_4 * contracted_part + 0.5 * gty_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Ttz(p.I) += d_t_Phi * d_z_Phi - 0.5 * gtz_4 * contracted_part + 0.5 * gtz_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Txx(p.I) += d_x_Phi * d_x_Phi - 0.5 * gxx_4 * contracted_part + 0.5 * gxx_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Txy(p.I) += d_x_Phi * d_y_Phi - 0.5 * gxy_4 * contracted_part + 0.5 * gxy_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Txz(p.I) += d_x_Phi * d_z_Phi - 0.5 * gxz_4 * contracted_part + 0.5 * gxz_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Tyy(p.I) += d_y_Phi * d_y_Phi - 0.5 * gyy_4 * contracted_part + 0.5 * gyy_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Tyz(p.I) += d_y_Phi * d_z_Phi - 0.5 * gyz_4 * contracted_part + 0.5 * gyz_4 * field_mass * field_mass * PhiL * PhiL;
    gf_Tzz(p.I) += d_z_Phi * d_z_Phi - 0.5 * gzz_4 * contracted_part + 0.5 * gzz_4 * field_mass * field_mass * PhiL * PhiL;
  };

  loop_int<0, 0, 0>(cctkGH, Tmunu_lambda);
}
} // namespace KleinGordonX
