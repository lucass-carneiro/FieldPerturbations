/*
 *  KleinGordon - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordon.
 *
 *  KleinGordon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Evolve.c
 * Implements the evolution equations for the ADM-decomposed scalar wave
 * equation as presented in Eqs. (A3c) and (A3d) of
 * https://arxiv.org/pdf/1709.06118.pdf.
 * Tensorial quantities produced in Wolfram Mathematica.
 * 4th order finite differencing
 */

/*************************
 * This thorn's includes *
 *************************/
#include "Derivatives.h"
#include "KleinGordon.h"

void KleinGordon_RHS_4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Ghost zone indexes */
  const CCTK_INT gx = cctk_nghostzones[0];
  const CCTK_INT gy = cctk_nghostzones[1];
  const CCTK_INT gz = cctk_nghostzones[2];

  /* Quantities required for the derivative macros to work */
  DECLARE_DERIVATIVE_FACTORS_4;

  /* Loop indexes */
  CCTK_INT i = 0, j = 0, k = 0, ijk = 0;

  /* Local ADM variables */
  CCTK_REAL alpL = 0.0;

  CCTK_REAL betaxL = 0.0;
  CCTK_REAL betayL = 0.0;
  CCTK_REAL betazL = 0.0;

  CCTK_REAL gxxL = 0.0;
  CCTK_REAL gxyL = 0.0;
  CCTK_REAL gxzL = 0.0;
  CCTK_REAL gyyL = 0.0;
  CCTK_REAL gyzL = 0.0;
  CCTK_REAL gzzL = 0.0;

  CCTK_REAL kxxL = 0.0;
  CCTK_REAL kxyL = 0.0;
  CCTK_REAL kxzL = 0.0;
  CCTK_REAL kyyL = 0.0;
  CCTK_REAL kyzL = 0.0;
  CCTK_REAL kzzL = 0.0;

  /* Local grid functions */
  CCTK_REAL PhiL = 0.0;
  CCTK_REAL K_PhiL = 0.0;

  /* Local Trace of extrinsic curvature */
  CCTK_REAL KTraceL = 0.0;

  /* Local Inverse metric */
  CCTK_REAL igxxL = 0.0;
  CCTK_REAL igxyL = 0.0;
  CCTK_REAL igxzL = 0.0;
  CCTK_REAL igyyL = 0.0;
  CCTK_REAL igyzL = 0.0;
  CCTK_REAL igzzL = 0.0;

  /* Local metric determinant */
  CCTK_REAL gdetL = 0.0;

  /* Parts of K_Phi right hand side */
  CCTK_REAL K_Phi_rhs_p1 = 0.0;
  CCTK_REAL K_Phi_rhs_p2 = 0.0;
  CCTK_REAL K_Phi_rhs_p3 = 0.0;
  CCTK_REAL K_Phi_rhs_p4 = 0.0;
  CCTK_REAL K_Phi_rhs_p5 = 0.0;

  /* Derivatives of Phi */
  CCTK_REAL d_x_Phi = 0.0;
  CCTK_REAL d_y_Phi = 0.0;
  CCTK_REAL d_z_Phi = 0.0;

  CCTK_REAL d_xx_Phi = 0.0;
  CCTK_REAL d_xy_Phi = 0.0;
  CCTK_REAL d_xz_Phi = 0.0;

  CCTK_REAL d_yy_Phi = 0.0;
  CCTK_REAL d_yz_Phi = 0.0;

  CCTK_REAL d_zz_Phi = 0.0;

  /* Derivatives of the metric */
  CCTK_REAL d_x_gxx = 0.0;
  CCTK_REAL d_y_gxx = 0.0;
  CCTK_REAL d_z_gxx = 0.0;

  CCTK_REAL d_x_gxy = 0.0;
  CCTK_REAL d_y_gxy = 0.0;
  CCTK_REAL d_z_gxy = 0.0;

  CCTK_REAL d_x_gxz = 0.0;
  CCTK_REAL d_y_gxz = 0.0;
  CCTK_REAL d_z_gxz = 0.0;

  CCTK_REAL d_x_gyy = 0.0;
  CCTK_REAL d_y_gyy = 0.0;
  CCTK_REAL d_z_gyy = 0.0;

  CCTK_REAL d_x_gyz = 0.0;
  CCTK_REAL d_y_gyz = 0.0;
  CCTK_REAL d_z_gyz = 0.0;

  CCTK_REAL d_x_gzz = 0.0;
  CCTK_REAL d_y_gzz = 0.0;
  CCTK_REAL d_z_gzz = 0.0;

  /* Derivatives of Alpha */
  CCTK_REAL d_x_alp = 0.0;
  CCTK_REAL d_y_alp = 0.0;
  CCTK_REAL d_z_alp = 0.0;

  /* Derivatives of K_Phi */
  CCTK_REAL d_x_K_Phi = 0.0;
  CCTK_REAL d_y_K_Phi = 0.0;
  CCTK_REAL d_z_K_Phi = 0.0;

  /* Christoffell symbols */
  CCTK_REAL Gamma_xxx = 0.0;
  CCTK_REAL Gamma_xxy = 0.0;
  CCTK_REAL Gamma_xxz = 0.0;
  CCTK_REAL Gamma_xyy = 0.0;
  CCTK_REAL Gamma_xyz = 0.0;
  CCTK_REAL Gamma_xzz = 0.0;

  CCTK_REAL Gamma_yxx = 0.0;
  CCTK_REAL Gamma_yxy = 0.0;
  CCTK_REAL Gamma_yxz = 0.0;
  CCTK_REAL Gamma_yyy = 0.0;
  CCTK_REAL Gamma_yyz = 0.0;
  CCTK_REAL Gamma_yzz = 0.0;

  CCTK_REAL Gamma_zxx = 0.0;
  CCTK_REAL Gamma_zxy = 0.0;
  CCTK_REAL Gamma_zxz = 0.0;
  CCTK_REAL Gamma_zyy = 0.0;
  CCTK_REAL Gamma_zyz = 0.0;
  CCTK_REAL Gamma_zzz = 0.0;

#pragma omp parallel for
  for (k = gz; k < cctk_lsh[2] - gz; k++) {
    for (j = gy; j < cctk_lsh[1] - gy; j++) {
      for (i = gx; i < cctk_lsh[0] - gx; i++) {
        ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        /* Assing ADM local variables */
        alpL = alp[ijk];

        betaxL = betax[ijk];
        betayL = betay[ijk];
        betazL = betaz[ijk];

        gxxL = gxx[ijk];
        gxyL = gxy[ijk];
        gxzL = gxz[ijk];
        gyyL = gyy[ijk];
        gyzL = gyz[ijk];
        gzzL = gzz[ijk];

        kxxL = kxx[ijk];
        kxyL = kxy[ijk];
        kxzL = kxz[ijk];
        kyyL = kyy[ijk];
        kyzL = kyz[ijk];
        kzzL = kzz[ijk];

        /* Assing wave eq. local variables */
        PhiL = Phi[ijk];
        K_PhiL = K_Phi[ijk];

        /* Computing the inverse metric */
        gdetL = -(gxzL * gxzL * gyyL) + 2 * gxyL * gxzL * gyzL -
                gxxL * gyzL * gyzL - gxyL * gxyL * gzzL + gxxL * gyyL * gzzL;
        igxxL = (-gyzL * gyzL + gyyL * gzzL) / gdetL;
        igxyL = (gxzL * gyzL - gxyL * gzzL) / gdetL;
        igxzL = (-(gxzL * gyyL) + gxyL * gyzL) / gdetL;
        igyyL = (-gxzL * gxzL + gxxL * gzzL) / gdetL;
        igyzL = (gxyL * gxzL - gxxL * gyzL) / gdetL;
        igzzL = (-gxyL * gxyL + gxxL * gyyL) / gdetL;

        /* Computing the trace of extrinsic curvature */
        KTraceL = igxxL * kxxL + igyyL * kyyL + igzzL * kzzL +
                  2 * igxyL * kxyL + 2 * igxzL * kxzL + 2 * igyzL * kyzL;

        /* Derivatives of Phi */
        d_x_Phi = D4x(Phi);
        d_y_Phi = D4y(Phi);
        d_z_Phi = D4z(Phi);

        d_xx_Phi = D4xx(Phi);
        d_xy_Phi = D4xy(Phi);
        d_xz_Phi = D4xz(Phi);

        d_yy_Phi = D4yy(Phi);
        d_yz_Phi = D4yz(Phi);

        d_zz_Phi = D4zz(Phi);

        /* Derivatives of the metric */
        d_x_gxx = D4x(gxx);
        d_y_gxx = D4y(gxx);
        d_z_gxx = D4z(gxx);

        d_x_gxy = D4x(gxy);
        d_y_gxy = D4y(gxy);
        d_z_gxy = D4z(gxy);

        d_x_gxz = D4x(gxz);
        d_y_gxz = D4y(gxz);
        d_z_gxz = D4z(gxz);

        d_x_gyy = D4x(gyy);
        d_y_gyy = D4y(gyy);
        d_z_gyy = D4z(gyy);

        d_x_gyz = D4x(gyz);
        d_y_gyz = D4y(gyz);
        d_z_gyz = D4z(gyz);

        d_x_gzz = D4x(gzz);
        d_y_gzz = D4y(gzz);
        d_z_gzz = D4z(gzz);

        /* Derivatives of Alpha */
        d_x_alp = D4x(alp);
        d_y_alp = D4y(alp);
        d_z_alp = D4z(alp);

        /* Derivatives of K_Phi */
        d_x_K_Phi = D4x(K_Phi);
        d_y_K_Phi = D4y(K_Phi);
        d_z_K_Phi = D4z(K_Phi);

        /* Christoffell symbols */
        Gamma_xxx = 0.5 * (igxxL * d_x_gxx - igxyL * d_y_gxx - igxzL * d_z_gxx +
                           2 * igxyL * d_x_gxy + 2 * igxzL * d_x_gxz);
        Gamma_xxy = 0.5 * (igxxL * d_y_gxx + igxyL * d_x_gyy +
                           igxzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
        Gamma_xxz =
            0.5 * (igxxL * d_z_gxx + igxyL * (d_z_gxy - d_y_gxz + d_x_gyz) +
                   igxzL * d_x_gzz);
        Gamma_xyy =
            0.5 * (2 * igxxL * d_y_gxy - igxxL * d_x_gyy + igxyL * d_y_gyy -
                   igxzL * d_z_gyy + 2 * igxzL * d_y_gyz);
        Gamma_xyz =
            0.5 * (igxyL * d_z_gyy + igxxL * (d_z_gxy + d_y_gxz - d_x_gyz) +
                   igxzL * d_y_gzz);
        Gamma_xzz = 0.5 * (2 * igxxL * d_z_gxz + 2 * igxyL * d_z_gyz -
                           igxxL * d_x_gzz - igxyL * d_y_gzz + igxzL * d_z_gzz);

        Gamma_yxx = 0.5 * (igxyL * d_x_gxx - igyyL * d_y_gxx - igyzL * d_z_gxx +
                           2 * igyyL * d_x_gxy + 2 * igyzL * d_x_gxz);
        Gamma_yxy = 0.5 * (igxyL * d_y_gxx + igyyL * d_x_gyy +
                           igyzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
        Gamma_yxz =
            0.5 * (igxyL * d_z_gxx + igyyL * (d_z_gxy - d_y_gxz + d_x_gyz) +
                   igyzL * d_x_gzz);
        Gamma_yyy =
            0.5 * (2 * igxyL * d_y_gxy - igxyL * d_x_gyy + igyyL * d_y_gyy -
                   igyzL * d_z_gyy + 2 * igyzL * d_y_gyz);
        Gamma_yyz =
            0.5 * (igyyL * d_z_gyy + igxyL * (d_z_gxy + d_y_gxz - d_x_gyz) +
                   igyzL * d_y_gzz);
        Gamma_yzz = 0.5 * (2 * igxyL * d_z_gxz + 2 * igyyL * d_z_gyz -
                           igxyL * d_x_gzz - igyyL * d_y_gzz + igyzL * d_z_gzz);

        Gamma_zxx = 0.5 * (igxzL * d_x_gxx - igyzL * d_y_gxx - igzzL * d_z_gxx +
                           2 * igyzL * d_x_gxy + 2 * igzzL * d_x_gxz);
        Gamma_zxy = 0.5 * (igxzL * d_y_gxx + igyzL * d_x_gyy +
                           igzzL * (-d_z_gxy + d_y_gxz + d_x_gyz));
        Gamma_zxz =
            0.5 * (igxzL * d_z_gxx + igyzL * (d_z_gxy - d_y_gxz + d_x_gyz) +
                   igzzL * d_x_gzz);
        Gamma_zyy =
            0.5 * (2 * igxzL * d_y_gxy - igxzL * d_x_gyy + igyzL * d_y_gyy -
                   igzzL * d_z_gyy + 2 * igzzL * d_y_gyz);
        Gamma_zyz =
            0.5 * (igyzL * d_z_gyy + igxzL * (d_z_gxy + d_y_gxz - d_x_gyz) +
                   igzzL * d_y_gzz);
        Gamma_zzz = 0.5 * (2 * igxzL * d_z_gxz + 2 * igyzL * d_z_gyz -
                           igxzL * d_x_gzz - igyzL * d_y_gzz + igzzL * d_z_gzz);

        /* Part 1 of K_Phi_rhs */
        K_Phi_rhs_p1 = KTraceL * K_PhiL;

        /* Part 2 of K_Phi_rhs */
        K_Phi_rhs_p2 =
            igxxL * d_xx_Phi + igxyL * (d_xy_Phi + d_xy_Phi) +
            igyyL * d_yy_Phi + igxzL * (d_xz_Phi + d_xz_Phi) +
            igyzL * (d_yz_Phi + d_yz_Phi) + igzzL * d_zz_Phi -
            d_x_Phi * (igxxL * Gamma_xxx + 2 * igxyL * Gamma_xxy +
                       2 * igxzL * Gamma_xxz + igyyL * Gamma_xyy +
                       2 * igyzL * Gamma_xyz + igzzL * Gamma_xzz) -
            igxxL * d_y_Phi * Gamma_yxx - 2 * igxyL * d_y_Phi * Gamma_yxy -
            2 * igxzL * d_y_Phi * Gamma_yxz - igyyL * d_y_Phi * Gamma_yyy -
            2 * igyzL * d_y_Phi * Gamma_yyz - igzzL * d_y_Phi * Gamma_yzz -
            igxxL * d_z_Phi * Gamma_zxx - 2 * igxyL * d_z_Phi * Gamma_zxy -
            2 * igxzL * d_z_Phi * Gamma_zxz - igyyL * d_z_Phi * Gamma_zyy -
            2 * igyzL * d_z_Phi * Gamma_zyz - igzzL * d_z_Phi * Gamma_zzz;

        /* Part 3 of K_Phi_rhs */
        K_Phi_rhs_p3 = field_mass * field_mass * PhiL;

        /* Part 4 of K_Phi_rhs */
        K_Phi_rhs_p4 =
            d_x_alp * (igxxL * d_x_Phi + igxyL * d_y_Phi + igxzL * d_z_Phi) +
            d_y_alp * (igxyL * d_x_Phi + igyyL * d_y_Phi + igyzL * d_z_Phi) +
            d_z_alp * (igxzL * d_x_Phi + igyzL * d_y_Phi + igzzL * d_z_Phi);

        /* Part 5 of K_Phi_rhs */
        K_Phi_rhs_p5 =
            betaxL * d_x_K_Phi + betayL * d_y_K_Phi + betazL * d_z_K_Phi;

        /* Phi_rhs */
        Phi_rhs[ijk] = -2.0 * alpL * K_PhiL + betaxL * d_x_Phi +
		  betayL * d_y_Phi + betazL * d_z_Phi;

        /* K_Phi_rhs */
        K_Phi_rhs[ijk] =
            alpL * (K_Phi_rhs_p1 - 0.5 * K_Phi_rhs_p2 + 0.5 * K_Phi_rhs_p3) -
            0.5 * K_Phi_rhs_p4 + K_Phi_rhs_p5;
      }
    }
  }
}
