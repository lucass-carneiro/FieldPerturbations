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
 *  along with KleinGordon.  If not, see <https://www.gnu.org/licenses/>.
 *
 * CalcTmunu_8.c
 * Compute the energy momentum tensor of the wave equation.
 */

/*************************
 * This thorn's includes *
 *************************/
#include "Derivatives.h"
#include "KleinGordon.h"

void KleinGordon_CalcTmunu_8(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Quantities required for the derivative macros to work */
  DECLARE_FIRST_DERIVATIVE_FACTORS_8;

  /* Loop index */
  CCTK_INT ijk = 0;

  /* Local ADM variables */
  CCTK_REAL alpL = 0.0;

  CCTK_REAL betaxL = 0.0;
  CCTK_REAL betayL = 0.0;
  CCTK_REAL betazL = 0.0;

  CCTK_REAL hxxL = 0.0;
  CCTK_REAL hxyL = 0.0;
  CCTK_REAL hxzL = 0.0;
  CCTK_REAL hyyL = 0.0;
  CCTK_REAL hyzL = 0.0;
  CCTK_REAL hzzL = 0.0;

  /* Local grid functions */
  CCTK_REAL PhiL = 0.0;
  CCTK_REAL K_PhiL = 0.0;

  /* Local Inverse metric */
  CCTK_REAL hdetL = 0.0;
  CCTK_REAL ihxxL = 0.0;
  CCTK_REAL ihxyL = 0.0;
  CCTK_REAL ihxzL = 0.0;
  CCTK_REAL ihyyL = 0.0;
  CCTK_REAL ihyzL = 0.0;
  CCTK_REAL ihzzL = 0.0;

  /* Local covariant (lower) shift vector */
  CCTK_REAL ibetaxL = 0.0;
  CCTK_REAL ibetayL = 0.0;
  CCTK_REAL ibetazL = 0.0;

  /*
   *Local 4-metric (lower)
   * It's only necessary to compute g_tt since the other componets are already
   * computed
   */
  CCTK_REAL gttL = 0.0;

  /* Inverse 4-metric */
  CCTK_REAL igttL = 0.0;
  CCTK_REAL igtxL = 0.0;
  CCTK_REAL igtyL = 0.0;
  CCTK_REAL igtzL = 0.0;
  CCTK_REAL igxxL = 0.0;
  CCTK_REAL igxyL = 0.0;
  CCTK_REAL igxzL = 0.0;
  CCTK_REAL igyyL = 0.0;
  CCTK_REAL igyzL = 0.0;
  CCTK_REAL igzzL = 0.0;

  /* Derivatives of Phi */
  CCTK_REAL d_x_Phi = 0.0;
  CCTK_REAL d_y_Phi = 0.0;
  CCTK_REAL d_z_Phi = 0.0;
  CCTK_REAL d_t_Phi = 0.0;

  // The scalar quantity g^{ab} \nabla_{a} \phi \nabla_{b} \phi
  CCTK_REAL nabladot = 0.0;

  /* Coordinate transformation jacobians */
  CCTK_REAL J11L = 0;
  CCTK_REAL J12L = 0;
  CCTK_REAL J13L = 0;

  CCTK_REAL J21L = 0;
  CCTK_REAL J22L = 0;
  CCTK_REAL J23L = 0;

  CCTK_REAL J31L = 0;
  CCTK_REAL J32L = 0;
  CCTK_REAL J33L = 0;

#pragma omp parallel
  CCTK_LOOP3_INT(loop_Tmunu, cctkGH, i, j, k) {
    ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    /* Assing ADM local variables */
    alpL = alp[ijk];

    betaxL = betax[ijk];
    betayL = betay[ijk];
    betazL = betaz[ijk];

    hxxL = gxx[ijk];
    hxyL = gxy[ijk];
    hxzL = gxz[ijk];
    hyyL = gyy[ijk];
    hyzL = gyz[ijk];
    hzzL = gzz[ijk];

    /* Assing wave eq. local variables */
    PhiL = Phi[ijk];
    K_PhiL = K_Phi[ijk];

    /* Assign Jacobias */
    J11L = J11[ijk];
    J12L = J12[ijk];
    J13L = J13[ijk];

    J21L = J21[ijk];
    J22L = J22[ijk];
    J23L = J23[ijk];

    J31L = J31[ijk];
    J32L = J32[ijk];
    J33L = J33[ijk];

    /* Computing the inverse 3-metric */
    hdetL = -(hxzL * hxzL * hyyL) + 2 * hxyL * hxzL * hyzL - hxxL * hyzL * hyzL - hxyL * hxyL * hzzL
            + hxxL * hyyL * hzzL;
    ihxxL = (-hyzL * hyzL + hyyL * hzzL) / hdetL;
    ihxyL = (hxzL * hyzL - hxyL * hzzL) / hdetL;
    ihxzL = (-(hxzL * hyyL) + hxyL * hyzL) / hdetL;
    ihyyL = (-hxzL * hxzL + hxxL * hzzL) / hdetL;
    ihyzL = (hxyL * hxzL - hxxL * hyzL) / hdetL;
    ihzzL = (-hxyL * hxyL + hxxL * hyyL) / hdetL;

    /* Computing the covariant (lower) shift vector */
    ibetaxL = hxxL * betaxL + hxyL * betayL + hxzL * betazL;
    ibetayL = hxyL * betaxL + hyyL * betayL + hyzL * betazL;
    ibetazL = hxzL * betaxL + hyzL * betayL + hzzL * betazL;

    /* Reconstructing the 4-metric (lower). */
    gttL = -(alpL * alpL) + ibetaxL * betaxL + ibetayL * betayL + ibetazL * betazL;

    // inverse 4-metric (upper)
    igttL = -1.0 / (alpL * alpL);
    igtxL = -1.0 * igttL * betaxL;
    igtyL = -1.0 * igttL * betayL;
    igtzL = -1.0 * igttL * betazL;
    igxxL = ihxxL + igttL * betaxL * betaxL;
    igxyL = ihxyL + igttL * betaxL * betayL;
    igxzL = ihxzL + igttL * betaxL * betazL;
    igyyL = ihyyL + igttL * betayL * betayL;
    igyzL = ihyzL + igttL * betayL * betazL;
    igzzL = ihzzL + igttL * betazL * betazL;

    /* Derivatives of Phi */
    d_x_Phi = global_Dx(8, Phi);
    d_y_Phi = global_Dy(8, Phi);
    d_z_Phi = global_Dz(8, Phi);
    d_t_Phi = (betaxL * d_x_Phi + betayL * d_y_Phi + betazL * d_z_Phi) - 2.0 * alpL * K_PhiL;

    // The scalar quantity g^{ab} \nabla_{a} \phi \nabla_{b} \phi
    nabladot = (igttL * d_t_Phi * d_t_Phi) + 2.0 * (igtxL * d_t_Phi * d_x_Phi)
               + 2.0 * (igtyL * d_t_Phi * d_y_Phi) + 2.0 * (igtzL * d_t_Phi * d_z_Phi)
               + (igxxL * d_x_Phi * d_x_Phi) + 2.0 * (igxyL * d_x_Phi * d_y_Phi)
               + 2.0 * (igxzL * d_x_Phi * d_z_Phi) + (igyyL * d_y_Phi * d_y_Phi)
               + 2.0 * (igyzL * d_y_Phi * d_z_Phi) + (igzzL * d_z_Phi * d_z_Phi);

    eTtt[ijk] += (d_t_Phi * d_t_Phi)
                 + 0.5 * gttL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTtx[ijk] += (d_t_Phi * d_x_Phi)
                 + 0.5 * betaxL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTty[ijk] += (d_t_Phi * d_y_Phi)
                 + 0.5 * betayL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTtz[ijk] += (d_t_Phi * d_z_Phi)
                 + 0.5 * betazL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTxx[ijk] += (d_x_Phi * d_x_Phi)
                 + 0.5 * hxxL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTxy[ijk] += (d_x_Phi * d_y_Phi)
                 + 0.5 * hxyL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTxz[ijk] += (d_x_Phi * d_z_Phi)
                 + 0.5 * hxzL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTyy[ijk] += (d_y_Phi * d_y_Phi)
                 + 0.5 * hyyL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTyz[ijk] += (d_y_Phi * d_z_Phi)
                 + 0.5 * hyzL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
    eTzz[ijk] += (d_z_Phi * d_z_Phi)
                 + 0.5 * hzzL * ((field_mass * PhiL) * (field_mass * PhiL) - nabladot);
  }
  CCTK_ENDLOOP3_INT(loop_Tmunu);
}
