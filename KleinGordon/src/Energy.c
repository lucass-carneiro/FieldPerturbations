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
#include "KleinGordon.h"
#include "Derivatives.h"

void KleinGordon_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Quantities required for the derivative macros to work */
  DECLARE_FIRST_DERIVATIVE_FACTORS_4;

  /* Loop indexes */
  CCTK_INT i = 0, j = 0, k = 0, ijk = 0;

  /* Local 3-metric, inverse 3-metric and 3-metric determinant*/
  CCTK_REAL gxxL = 0.0;
  CCTK_REAL gxyL = 0.0;
  CCTK_REAL gxzL = 0.0;
  CCTK_REAL gyyL = 0.0;
  CCTK_REAL gyzL = 0.0;
  CCTK_REAL gzzL = 0.0;

  CCTK_REAL igxxL = 0.0;
  CCTK_REAL igxyL = 0.0;
  CCTK_REAL igxzL = 0.0;
  CCTK_REAL igyyL = 0.0;
  CCTK_REAL igyzL = 0.0;
  CCTK_REAL igzzL = 0.0;

  CCTK_REAL gdetL = 0.0;

  /* Derivatives of Phi */
  CCTK_REAL d_x_Phi = 0.0;
  CCTK_REAL d_y_Phi = 0.0;
  CCTK_REAL d_z_Phi = 0.0;

  /* Parts of the energy density */
  CCTK_REAL contracted_part = 0;

#pragma omp parallel for
  for (k = 0; k < cctk_lsh[2]; k++) {
    for (j = 0; j < cctk_lsh[1]; j++) {
      for (i = 0; i < cctk_lsh[0]; i++) {
        ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        /* The 3-metric */
        gxxL = gxx[ijk];
        gxyL = gxy[ijk];
        gxzL = gxz[ijk];
        gyyL = gyy[ijk];
        gyzL = gyz[ijk];
        gzzL = gzz[ijk];

        /* The 3-metric inverse */
        gdetL = -(gxzL * gxzL * gyyL) + 2 * gxyL * gxzL * gyzL -
                gxxL * gyzL * gyzL - gxyL * gxyL * gzzL + gxxL * gyyL * gzzL;
        igxxL = (-gyzL * gyzL + gyyL * gzzL) / gdetL;
        igxyL = (gxzL * gyzL - gxyL * gzzL) / gdetL;
        igxzL = (-(gxzL * gyyL) + gxyL * gyzL) / gdetL;
        igyyL = (-gxzL * gxzL + gxxL * gzzL) / gdetL;
        igyzL = (gxyL * gxzL - gxxL * gyzL) / gdetL;
        igzzL = (-gxyL * gxyL + gxxL * gyyL) / gdetL;

        /* The derivatives of Phi */
        d_x_Phi = D4x(Phi);
        d_y_Phi = D4y(Phi);
        d_z_Phi = D4z(Phi);

        /* Contracted part */
        contracted_part =
            igxxL * d_x_Phi * d_x_Phi + igxyL * d_x_Phi * d_y_Phi +
            igxzL * d_x_Phi * d_z_Phi + igyyL * d_y_Phi * d_y_Phi +
            igyzL * d_y_Phi * d_z_Phi + igzzL * d_z_Phi * d_z_Phi;

        epsilon[ijk] = 0.5 * (K_Phi[ijk] * K_Phi[ijk] + contracted_part);
      }
    }
  }
}
