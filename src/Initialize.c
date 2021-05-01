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
 *  Initialize.c
 *  Initialize grid variables.
 */

/*******************
 * Cactus includes *
 *******************/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/************************
 * C std. lib. includes *
 ************************/
#include <math.h>

/*****************
 * Local Defines *
 *****************/
#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

/**************
 * Prototypes *
 **************/
void ADMScalarWave_Initialize(CCTK_ARGUMENTS);

/**************************************************
 * ADMScalarWave_Initialize(CCTK_ARGUMENTS)       *
 *                                                *
 * This function provies the scalar fields with   *
 * initial data to begin it's time evolution.     *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void ADMScalarWave_Initialize(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Loop indexes */
  int i = 0;
  int j = 0;
  int k = 0;

  /* Determine which type of initial data to apply */
  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    CCTK_REAL xmx0 = 0.0;
    CCTK_REAL ymy0 = 0.0;
    CCTK_REAL zmz0 = 0.0;

    CCTK_REAL xmx02 = 0.0;
    CCTK_REAL ymy02 = 0.0;
    CCTK_REAL zmz02 = 0.0;

    CCTK_REAL R2 = 0.0;
    CCTK_REAL R = 0.0;

    CCTK_REAL R2_shifted = 0.0;
    CCTK_REAL R_shifted = 0.0;

    CCTK_REAL expo = 0.0;

    CCTK_REAL dipole = 0.0;
    CCTK_REAL quadrupole = 0.0;

    CCTK_REAL amplitude = 0.0;

    size_t ijk = 0;

    /* Loop over all points (ghostzones included) */
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {

          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          xmx0 = x[ijk] - gaussian_x0;
          ymy0 = y[ijk] - gaussian_y0;
          zmz0 = z[ijk] - gaussian_z0;

          xmx02 = SQR(xmx0);
          ymy02 = SQR(ymy0);
          zmz02 = SQR(zmz0);

          R2 = xmx02 + ymy02 + zmz02;
          R = sqrt(R2);

          R2_shifted = R2 + 1.0e-5;
          R_shifted = sqrt(R2_shifted);

          expo = exp(-0.5 * SQR((R - gaussian_R0) / gaussian_sigma));

          dipole = (xmx0 - zmz0) / R_shifted;
          quadrupole = (xmx02 - ymy02 + zmz02 + xmx0 * zmz0) / R2_shifted;

          amplitude =
              gaussian_c0 + gaussian_c1 * dipole + gaussian_c2 * quadrupole;

          Phi[ijk] = amplitude * expo;
          K_Phi[ijk] = 0.0;

          Phi_rhs[ijk] = 0.0;
          K_Phi_rhs[ijk] = 0.0;
        }
      }
    }
  }
}
