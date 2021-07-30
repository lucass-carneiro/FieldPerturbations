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
 * Energy.c
 * Calculate the wave equation's associated energy, as
 * presented in https://arxiv.org/abs/1005.2922v1.
 */

/*************************
 * This thorn's includes *
 *************************/
#include "Derivatives.h"
#include "KleinGordon.h"

void KleinGordon_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Quantities required for the derivative macros to work */
  DECLARE_FIRST_DERIVATIVE_FACTORS_4;

  CCTK_LOOP3_INT(loop_energy, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    /* The 3-metric */
    const CCTK_REAL gxxL = gxx[ijk];
    const CCTK_REAL gxyL = gxy[ijk];
    const CCTK_REAL gxzL = gxz[ijk];
    const CCTK_REAL gyyL = gyy[ijk];
    const CCTK_REAL gyzL = gyz[ijk];
    const CCTK_REAL gzzL = gzz[ijk];

    /* The 3-metric inverse */
    const CCTK_REAL gdetL = -(gxzL * gxzL * gyyL) + 2 * gxyL * gxzL * gyzL -
                            gxxL * gyzL * gyzL - gxyL * gxyL * gzzL +
                            gxxL * gyyL * gzzL;
    const CCTK_REAL igxxL = (-gyzL * gyzL + gyyL * gzzL) / gdetL;
    const CCTK_REAL igxyL = (gxzL * gyzL - gxyL * gzzL) / gdetL;
    const CCTK_REAL igxzL = (-(gxzL * gyyL) + gxyL * gyzL) / gdetL;
    const CCTK_REAL igyyL = (-gxzL * gxzL + gxxL * gzzL) / gdetL;
    const CCTK_REAL igyzL = (gxyL * gxzL - gxxL * gyzL) / gdetL;
    const CCTK_REAL igzzL = (-gxyL * gxyL + gxxL * gyyL) / gdetL;

    /* The derivatives of Phi */
    const CCTK_REAL d_x_Phi = D4x(Phi);
    const CCTK_REAL d_y_Phi = D4y(Phi);
    const CCTK_REAL d_z_Phi = D4z(Phi);

    /* Contracted part */
    const CCTK_REAL contracted_part =
        igxxL * d_x_Phi * d_x_Phi + igxyL * d_x_Phi * d_y_Phi +
        igxzL * d_x_Phi * d_z_Phi + igyyL * d_y_Phi * d_y_Phi +
        igyzL * d_y_Phi * d_z_Phi + igzzL * d_z_Phi * d_z_Phi;

    epsilon[ijk] = 0.5 * (K_Phi[ijk] * K_Phi[ijk] + contracted_part);
  }
  CCTK_ENDLOOP3_INT(loop_energy);
}
