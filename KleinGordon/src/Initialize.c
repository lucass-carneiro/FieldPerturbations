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
void KleinGordon_Initialize(CCTK_ARGUMENTS);
CCTK_REAL exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);
CCTK_REAL dt_exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);
CCTK_REAL multipolar_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL);

/**************************************************
 * KleinGordon_Initialize(CCTK_ARGUMENTS)       *
 *                                                *
 * This function provies the scalar fields with   *
 * initial data to begin it's time evolution.     *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void KleinGordon_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Loop indexes */
  CCTK_INT i = 0, j = 0, k = 0, ijk = 0;

  /* Determine which type of initial data to apply */
  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    /* Loop over all points (ghostzones included) */
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          Phi[ijk] = multipolar_gaussian(x[ijk], y[ijk], z[ijk]);
          K_Phi[ijk] = 0.0;
        }
      }
    }
  } else if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

    /* Loop over all points (ghostzones included) */
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          Phi[ijk] = exact_gaussian(0.0, x[ijk], y[ijk], z[ijk]);
          K_Phi[ijk] = dt_exact_gaussian(0.0, x[ijk], y[ijk], z[ijk]);
        }
      }
    }
  }
}

inline CCTK_REAL f(CCTK_REAL x, CCTK_REAL sigma) {
  return exp(-0.5 * (x / sigma) * (x / sigma));
}
inline CCTK_REAL fx(CCTK_REAL x, CCTK_REAL sigma) {
  return (-x * f(x, sigma)) / (sigma * sigma);
}

/**********************************************************************
 * exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                                    *
 * Computes the exect result of a gaussian fild in a Minkowski        *
 * background for a given time and position                           *
 *                                                                    *
 * Input: The 4-D point where the gaussian shoulde be computed        *
 *                                                                    *
 * Output: The result of the gaussian at the given point              *
 **********************************************************************/
CCTK_REAL exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL r = sqrt((x - gaussian_x0) * (x - gaussian_x0) +
                           (y - gaussian_y0) * (y - gaussian_y0) +
                           (z - gaussian_z0) * (z - gaussian_z0));

  // Use L'Hôpital's rule for small r
  if (r < 1.0e-8)
    return fx(r - t, gaussian_sigma) - fx(r + t, gaussian_sigma);
  else
    return (f(r - t, gaussian_sigma) - f(r + t, gaussian_sigma)) / r;
}

/*************************************************************************
 * dt_exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                                       *
 * Computes the time derivative exect result of a gaussian field in a    *
 * Minkowski background for a given time and position                    *
 *                                                                       *
 * Input: The 4-D point where the time derivative if the gaussian should *
 * be computed                                                           *
 *                                                                       *
 * Output: The result of the time derivative of the gaussian at the      *
 * given point                                                           *
 *************************************************************************/
CCTK_REAL dt_exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL r = sqrt((x - gaussian_x0) * (x - gaussian_x0) +
                           (y - gaussian_y0) * (y - gaussian_y0) +
                           (z - gaussian_z0) * (z - gaussian_z0));

  // Use L'Hôpital's rule for small r
  if (r < 1.0e-8) {
    return (2 * (-t + gaussian_sigma) * (t + gaussian_sigma)) /
           (exp(t * t / (2. * gaussian_sigma * gaussian_sigma)) *
            gaussian_sigma * gaussian_sigma * gaussian_sigma * gaussian_sigma);
  } else {
    return (r + exp((2 * r * t) / (gaussian_sigma * gaussian_sigma)) * (r - t) +
            t) /
           (exp(((r + t) * (r + t)) / (2 * gaussian_sigma * gaussian_sigma)) *
            r * gaussian_sigma * gaussian_sigma);
  }
}

/**********************************************************************
 * exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                                    *
 * Computes the multipolar gaussian of zilhao                         *
 *                                                                    *
 * Input: The 4-D point where the gaussian shoulde be computed        *
 *                                                                    *
 * Output: The result of the gaussian at the given point              *
 **********************************************************************/
CCTK_REAL multipolar_gaussian(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL xmx0 = x - gaussian_x0;
  CCTK_REAL ymy0 = y - gaussian_y0;
  CCTK_REAL zmz0 = z - gaussian_z0;

  CCTK_REAL xmx02 = xmx0 * xmx0;
  CCTK_REAL ymy02 = ymy0 * ymy0;
  CCTK_REAL zmz02 = zmz0 * zmz0;

  CCTK_REAL R2 = xmx02 + ymy02 + zmz02;
  CCTK_REAL R = sqrt(R2);

  CCTK_REAL R2_shifted = R2 + 1.0e-5;
  CCTK_REAL R_shifted = sqrt(R2_shifted);

  CCTK_REAL expo = exp(-0.5 * SQR((R - gaussian_R0) / gaussian_sigma));

  CCTK_REAL dipole = (xmx0 - zmz0) / R_shifted;
  CCTK_REAL quadrupole = (xmx02 - ymy02 + zmz02 + xmx0 * zmz0) / R2_shifted;

  return gaussian_c0 + gaussian_c1 * dipole + gaussian_c2 * quadrupole;
}
