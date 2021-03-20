/*
 *  ADMScalarWave - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of ADMScalarWave.
 *
 *  ADMScalarWave is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADMScalarWave is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * ExcisionMask.c
 * Implements simple excision mask functions
 */

/*******************
 * Cactus includes *
 *******************/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/*****************
 * Libc includes *
 *****************/
#include <math.h>

/**************
 * Prototypes *
 **************/
void ADMScalarWave_ExcisionMask(CCTK_ARGUMENTS);
inline CCTK_REAL exp_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh);
inline CCTK_REAL cos_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh);
inline CCTK_REAL hermite_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh);
inline CCTK_REAL cubic_step(CCTK_REAL r, CCTK_REAL rh);

/**********************************************
 * ADMScalarWave_ExcisionMask(CCTK_ARGUMENTS) *
 *                                            *
 * This function implements simple excision   *
 * masks M(r) such that                       *
 *                                            *
 *        1 if r >= rh                        *
 * M(r) = 0 if r <= r0                        *
 *        g(r) if r0 < r < rh                 *
 *                                            *
 * where the g(r) profile is choosen by the   *
 * user.                                      *
 *                                            *
 * This mask function is multiplied by all    *
 * evolved functions AFTER the RHS is         *
 * calculated and boundary conditions are     *
 * applied.                                   *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void ADMScalarWave_ExcisionMask(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rL = 0.0;
  CCTK_REAL PhiL = 0.0;
  CCTK_REAL K_PhiL = 0.0;

  int i = 0, j = 0, k = 0, ijk = 0;

  if (CCTK_EQUALS(mask_type, "exp_smooth_step")) {
/* Loop over all points (ghostzones included) */
#pragma omp parallel for
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {

          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          rL = r[ijk];
          PhiL = Phi[ijk];
          K_PhiL = K_Phi[ijk];

          Phi[ijk] =
              exp_smooth_step(rL, excision_min_radius, excision_max_radius) *
              PhiL;
          K_Phi[ijk] =
              exp_smooth_step(rL, excision_min_radius, excision_max_radius) *
              K_PhiL;
        }
      }
    }
  } else if (CCTK_EQUALS(mask_type, "cos_smooth_step")) {
/* Loop over all points (ghostzones included) */
#pragma omp parallel for
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {

          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          rL = r[ijk];
          PhiL = Phi[ijk];
          K_PhiL = K_Phi[ijk];

          Phi[ijk] =
              cos_smooth_step(rL, excision_min_radius, excision_max_radius) *
              PhiL;
          K_Phi[ijk] =
              cos_smooth_step(rL, excision_min_radius, excision_max_radius) *
              K_PhiL;
        }
      }
    }
  } else if (CCTK_EQUALS(mask_type, "cubic_step")) {
/* Loop over all points (ghostzones included) */
#pragma omp parallel for
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {

          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          rL = r[ijk];
          PhiL = Phi[ijk];
          K_PhiL = K_Phi[ijk];

          Phi[ijk] = cubic_step(rL, excision_max_radius) * PhiL;
          K_Phi[ijk] = cubic_step(rL, excision_max_radius) * K_PhiL;
        }
      }
    }
  } else if (CCTK_EQUALS(mask_type, "hermite_smooth_step")) {
/* Loop over all points (ghostzones included) */
#pragma omp parallel for
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {

          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          rL = r[ijk];
          PhiL = Phi[ijk];
          K_PhiL = K_Phi[ijk];

          Phi[ijk] = hermite_smooth_step(rL, excision_min_radius,
                                         excision_max_radius) *
                     PhiL;
          K_Phi[ijk] = hermite_smooth_step(rL, excision_min_radius,
                                           excision_max_radius) *
                       K_PhiL;
        }
      }
    }
  }
}

/************************************************************
 * exp_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL r0) *
 *                                                          *
 * Implement the smooth step function as described in       *
 * https://en.wikipedia.org/wiki/ \                         *
 * Non-analytic_smooth_function#Smooth_transition_functions *
 *                                                          *
 * Input:                                                   *
 * CCTK_REAL r: The point of evaluation                     *
 * CCTK_REAL r0: The initial excision radius                *
 * CCTK_REAL rh: The final excision radius                  *
 * from interface.ccl                                       *
 *                                                          *
 * Output: Nothing                                          *
 ************************************************************/
inline CCTK_REAL exp_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh) {
  const CCTK_REAL y = (r - r0) / (rh - r0);
  const CCTK_REAL omy = 1.0 - y;

  const CCTK_REAL fy = (y > 0) ? exp(-1.0 / y) : 0.0;
  const CCTK_REAL fomy = (omy > 0) ? exp(-1.0 / omy) : 0.0;

  return fy / (fy + fomy);
}

/************************************************************
 * cos_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL r0) *
 *                                                          *
 * Implement the smooth step function as described using    *
 * the cossine function as the transition profile           *
 *                                                          *
 * Input:                                                   *
 * CCTK_REAL r: The point of evaluation                     *
 * CCTK_REAL r0: The initial excision radius                *
 * CCTK_REAL rh: The final excision radius                  *
 * from interface.ccl                                       *
 *                                                          *
 * Output: Nothing                                          *
 ************************************************************/
inline CCTK_REAL cos_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh) {
  if (r < r0)
    return 0.0;
  else if (r > rh)
    return 1.0;
  else
    return 0.5 + 0.5 * cos((M_PI / (rh - r0)) * (r - r0) + M_PI);
}

/*********************************************
 * cubic_step(CCTK_REAL r, CCTK_REAL r0)     *
 *                                           *
 * Implement a non smooth step function      *
 * using a cubic polynomial as the           *
 * transition profile                        *
 *                                           *
 * Input:                                    *
 * CCTK_REAL r: The point of evaluation      *
 * CCTK_REAL r0: The initial excision radius *
 * CCTK_REAL rh: The final excision radius   *
 * from interface.ccl                        *
 *                                           *
 * Output: Nothing                           *
 *********************************************/
inline CCTK_REAL cubic_step(CCTK_REAL r, CCTK_REAL rh) {
  return (r < rh) ? (r / rh) * (r / rh) * (r / rh) : 1.0;
}

/****************************************************************
 * hermite_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL r0) *
 *                                                              *
 * Implement the smooth step function as described using        *
 * a cubic hermite polynomial as the transition profile         *
 *                                                              *
 * Input:                                                       *
 * CCTK_REAL r: The point of evaluation                         *
 * CCTK_REAL r0: The initial excision radius                    *
 * CCTK_REAL rh: The final excision radius                      *
 * from interface.ccl                                           *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
inline CCTK_REAL hermite_smooth_step(CCTK_REAL r, CCTK_REAL r0, CCTK_REAL rh) {
  if (r > rh)
    return 1.0;
  else if (r < rh)
    return 0.0;
  else {
    CCTK_REAL y = (r - r0) / (rh - r0);
    y = y * y * (3.0 - 2.0 * y);
    return y;
  }
}
