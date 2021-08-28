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
 * RHS_6.cpp
 * Compute the RHS of the wave equation using 6th order accurate stencils
 */

#include "Derivative.hpp"
#include "KleinGordonX.hpp"

namespace KleinGordonX {

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_int;
using Loop::PointDesc;

extern "C" void KleinGordonX_RHS_6(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_RHS_6;
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

  const GF3D2<const CCTK_REAL> gf_kxx(layout, kxx);
  const GF3D2<const CCTK_REAL> gf_kxy(layout, kxy);
  const GF3D2<const CCTK_REAL> gf_kxz(layout, kxz);
  const GF3D2<const CCTK_REAL> gf_kyy(layout, kyy);
  const GF3D2<const CCTK_REAL> gf_kyz(layout, kyz);
  const GF3D2<const CCTK_REAL> gf_kzz(layout, kzz);

  // Right hand side ---------------------------------------------
  const GF3D2<CCTK_REAL> gf_Phi_rhs(layout, Phi_rhs);
  const GF3D2<CCTK_REAL> gf_K_Phi_rhs(layout, K_Phi_rhs);

  // Derivative objects -----------------------------------------
  const Derivative<Stencil::c7> d_Phi(gf_Phi);
  const Derivative<Stencil::c7> d_K_Phi(gf_K_Phi);
  const Derivative<Stencil::c7> d_alp(gf_alp);
  const Derivative<Stencil::c7> d_gxx(gf_gxx);
  const Derivative<Stencil::c7> d_gxy(gf_gxy);
  const Derivative<Stencil::c7> d_gxz(gf_gxz);
  const Derivative<Stencil::c7> d_gyy(gf_gyy);
  const Derivative<Stencil::c7> d_gyz(gf_gyz);
  const Derivative<Stencil::c7> d_gzz(gf_gzz);

  auto rhs_lambda = [&](const PointDesc &p) {
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

    const CCTK_REAL kxxL = gf_kxx(p.I);
    const CCTK_REAL kxyL = gf_kxy(p.I);
    const CCTK_REAL kxzL = gf_kxz(p.I);
    const CCTK_REAL kyyL = gf_kyy(p.I);
    const CCTK_REAL kyzL = gf_kyz(p.I);
    const CCTK_REAL kzzL = gf_kzz(p.I);

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

    // Computing the trace of extrinsic curvature
    const CCTK_REAL KTraceL = igxxL * kxxL + igyyL * kyyL + igzzL * kzzL +
                              2 * igxyL * kxyL + 2 * igxzL * kxzL +
                              2 * igyzL * kyzL;

    // Derivatives of Phi
    auto grad_Phi = d_Phi.grad(p);
    auto hess_Phi = d_Phi.hess(p);

    const CCTK_REAL d_x_Phi = grad_Phi[0];
    const CCTK_REAL d_y_Phi = grad_Phi[1];
    const CCTK_REAL d_z_Phi = grad_Phi[2];

    const CCTK_REAL d_xx_Phi = hess_Phi[0];
    const CCTK_REAL d_yy_Phi = hess_Phi[1];
    const CCTK_REAL d_zz_Phi = hess_Phi[2];

    const CCTK_REAL d_xy_Phi = hess_Phi[3];
    const CCTK_REAL d_xz_Phi = hess_Phi[4];
    const CCTK_REAL d_yz_Phi = hess_Phi[5];

    // Derivatives of the metric
    auto grad_gxx = d_gxx.grad(p);
    auto grad_gxy = d_gxy.grad(p);
    auto grad_gxz = d_gxz.grad(p);
    auto grad_gyy = d_gyy.grad(p);
    auto grad_gyz = d_gyz.grad(p);
    auto grad_gzz = d_gzz.grad(p);

    const CCTK_REAL d_x_gxx = grad_gxx[0];
    const CCTK_REAL d_y_gxx = grad_gxx[1];
    const CCTK_REAL d_z_gxx = grad_gxx[2];

    const CCTK_REAL d_x_gxy = grad_gxy[0];
    const CCTK_REAL d_y_gxy = grad_gxy[1];
    const CCTK_REAL d_z_gxy = grad_gxy[2];

    const CCTK_REAL d_x_gxz = grad_gxz[0];
    const CCTK_REAL d_y_gxz = grad_gxz[1];
    const CCTK_REAL d_z_gxz = grad_gxz[2];

    const CCTK_REAL d_x_gyy = grad_gyy[0];
    const CCTK_REAL d_y_gyy = grad_gyy[1];
    const CCTK_REAL d_z_gyy = grad_gyy[2];

    const CCTK_REAL d_x_gyz = grad_gyz[0];
    const CCTK_REAL d_y_gyz = grad_gyz[1];
    const CCTK_REAL d_z_gyz = grad_gyz[2];

    const CCTK_REAL d_x_gzz = grad_gzz[0];
    const CCTK_REAL d_y_gzz = grad_gzz[1];
    const CCTK_REAL d_z_gzz = grad_gzz[2];

    // Derivatives of Alpha
    auto grad_alp = d_alp.grad(p);

    const CCTK_REAL d_x_alp = grad_alp[0];
    const CCTK_REAL d_y_alp = grad_alp[1];
    const CCTK_REAL d_z_alp = grad_alp[2];

    // Derivatives of K_Phi
    auto grad_K_Phi = d_K_Phi.grad(p);

    const CCTK_REAL d_x_K_Phi = grad_K_Phi[0];
    const CCTK_REAL d_y_K_Phi = grad_K_Phi[1];
    const CCTK_REAL d_z_K_Phi = grad_K_Phi[2];

    // Christoffell symbols */
    const CCTK_REAL Gamma_xxx =
        0.5 * (igxxL * d_x_gxx - igxyL * d_y_gxx - igxzL * d_z_gxx +
               2 * igxyL * d_x_gxy + 2 * igxzL * d_x_gxz);
    const CCTK_REAL Gamma_xxy = 0.5 * (igxxL * d_y_gxx + igxyL * d_x_gyy +
                                       igxzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
    const CCTK_REAL Gamma_xxz =
        0.5 * (igxxL * d_z_gxx + igxyL * (d_z_gxy - d_y_gxz + d_x_gyz) +
               igxzL * d_x_gzz);
    const CCTK_REAL Gamma_xyy =
        0.5 * (2 * igxxL * d_y_gxy - igxxL * d_x_gyy + igxyL * d_y_gyy -
               igxzL * d_z_gyy + 2 * igxzL * d_y_gyz);
    const CCTK_REAL Gamma_xyz =
        0.5 * (igxyL * d_z_gyy + igxxL * (d_z_gxy + d_y_gxz - d_x_gyz) +
               igxzL * d_y_gzz);
    const CCTK_REAL Gamma_xzz =
        0.5 * (2 * igxxL * d_z_gxz + 2 * igxyL * d_z_gyz - igxxL * d_x_gzz -
               igxyL * d_y_gzz + igxzL * d_z_gzz);

    const CCTK_REAL Gamma_yxx =
        0.5 * (igxyL * d_x_gxx - igyyL * d_y_gxx - igyzL * d_z_gxx +
               2 * igyyL * d_x_gxy + 2 * igyzL * d_x_gxz);
    const CCTK_REAL Gamma_yxy = 0.5 * (igxyL * d_y_gxx + igyyL * d_x_gyy +
                                       igyzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
    const CCTK_REAL Gamma_yxz =
        0.5 * (igxyL * d_z_gxx + igyyL * (d_z_gxy - d_y_gxz + d_x_gyz) +
               igyzL * d_x_gzz);
    const CCTK_REAL Gamma_yyy =
        0.5 * (2 * igxyL * d_y_gxy - igxyL * d_x_gyy + igyyL * d_y_gyy -
               igyzL * d_z_gyy + 2 * igyzL * d_y_gyz);
    const CCTK_REAL Gamma_yyz =
        0.5 * (igyyL * d_z_gyy + igxyL * (d_z_gxy + d_y_gxz - d_x_gyz) +
               igyzL * d_y_gzz);
    const CCTK_REAL Gamma_yzz =
        0.5 * (2 * igxyL * d_z_gxz + 2 * igyyL * d_z_gyz - igxyL * d_x_gzz -
               igyyL * d_y_gzz + igyzL * d_z_gzz);

    const CCTK_REAL Gamma_zxx =
        0.5 * (igxzL * d_x_gxx - igyzL * d_y_gxx - igzzL * d_z_gxx +
               2 * igyzL * d_x_gxy + 2 * igzzL * d_x_gxz);
    const CCTK_REAL Gamma_zxy = 0.5 * (igxzL * d_y_gxx + igyzL * d_x_gyy +
                                       igzzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
    const CCTK_REAL Gamma_zxz =
        0.5 * (igxzL * d_z_gxx + igyzL * (d_z_gxy - d_y_gxz + d_x_gyz) +
               igzzL * d_x_gzz);
    const CCTK_REAL Gamma_zyy =
        0.5 * (2 * igxzL * d_y_gxy - igxzL * d_x_gyy + igyzL * d_y_gyy -
               igzzL * d_z_gyy + 2 * igzzL * d_y_gyz);
    const CCTK_REAL Gamma_zyz =
        0.5 * (igyzL * d_z_gyy + igxzL * (d_z_gxy + d_y_gxz - d_x_gyz) +
               igzzL * d_y_gzz);
    const CCTK_REAL Gamma_zzz =
        0.5 * (2 * igxzL * d_z_gxz + 2 * igyzL * d_z_gyz - igxzL * d_x_gzz -
               igyzL * d_y_gzz + igzzL * d_z_gzz);

    // Part 1 of K_Phi_rhs
    const CCTK_REAL K_Phi_rhs_p1 = KTraceL * K_PhiL;

    // Part 2 of const CCTK_REAL K_Phi_rhs
    const CCTK_REAL K_Phi_rhs_p2 =
        igxxL * d_xx_Phi + igxyL * (d_xy_Phi + d_xy_Phi) + igyyL * d_yy_Phi +
        igxzL * (d_xz_Phi + d_xz_Phi) + igyzL * (d_yz_Phi + d_yz_Phi) +
        igzzL * d_zz_Phi -
        d_x_Phi *
            (igxxL * Gamma_xxx + 2 * igxyL * Gamma_xxy + 2 * igxzL * Gamma_xxz +
             igyyL * Gamma_xyy + 2 * igyzL * Gamma_xyz + igzzL * Gamma_xzz) -
        igxxL * d_y_Phi * Gamma_yxx - 2 * igxyL * d_y_Phi * Gamma_yxy -
        2 * igxzL * d_y_Phi * Gamma_yxz - igyyL * d_y_Phi * Gamma_yyy -
        2 * igyzL * d_y_Phi * Gamma_yyz - igzzL * d_y_Phi * Gamma_yzz -
        igxxL * d_z_Phi * Gamma_zxx - 2 * igxyL * d_z_Phi * Gamma_zxy -
        2 * igxzL * d_z_Phi * Gamma_zxz - igyyL * d_z_Phi * Gamma_zyy -
        2 * igyzL * d_z_Phi * Gamma_zyz - igzzL * d_z_Phi * Gamma_zzz;

    // Part 3 of const CCTK_REAL K_Phi_rhs
    const CCTK_REAL K_Phi_rhs_p3 = field_mass * field_mass * PhiL;

    // Part 4 of const CCTK_REAL K_Phi_rhs */
    const CCTK_REAL K_Phi_rhs_p4 =
        d_x_alp * (igxxL * d_x_Phi + igxyL * d_y_Phi + igxzL * d_z_Phi) +
        d_y_alp * (igxyL * d_x_Phi + igyyL * d_y_Phi + igyzL * d_z_Phi) +
        d_z_alp * (igxzL * d_x_Phi + igyzL * d_y_Phi + igzzL * d_z_Phi);

    // Part 5 of const CCTK_REAL K_Phi_rhs */
    const CCTK_REAL K_Phi_rhs_p5 =
        betaxL * d_x_K_Phi + betayL * d_y_K_Phi + betazL * d_z_K_Phi;

    // Phi_rhs
    gf_Phi_rhs(p.I) = -2.0 * alpL * K_PhiL + betaxL * d_x_Phi +
                      betayL * d_y_Phi + betazL * d_z_Phi;
    // K_Phi_rhs
    gf_K_Phi_rhs(p.I) =
        alpL * (K_Phi_rhs_p1 - 0.5 * K_Phi_rhs_p2 + 0.5 * K_Phi_rhs_p3) -
        0.5 * K_Phi_rhs_p4 + K_Phi_rhs_p5;
  };

  loop_int<0, 0, 0>(cctkGH, rhs_lambda);
}
} // namespace KleinGordonX
