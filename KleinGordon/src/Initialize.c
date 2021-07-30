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

/*************************
 * This thorn's includes *
 *************************/
#include "Derivatives.h"
#include "KleinGordon.h"

/**************************
 * C std. lib. includes   *
 * and external libraries *
 **************************/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

inline CCTK_REAL f(CCTK_REAL x, CCTK_REAL sigma) {
  return exp(-0.5 * (x / sigma) * (x / sigma));
}
inline CCTK_REAL fx(CCTK_REAL x, CCTK_REAL sigma) {
  return (-x * f(x, sigma)) / (sigma * sigma);
}

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

/* This function allocates a buffer for storing the spherical harmonics
 * IT'S ASSUMED THAT THE BUFFER OWNERSHIP IS MOVED ON ASSIGNMENT,
 * thus the final user is responsible for freeing the memory
 * If memory allocation fails, the function halts Cactus
 */
CCTK_REAL *Ylm_buffer(CCTK_INT lmax) {

  CCTK_REAL *buffer = malloc(gsl_sf_legendre_array_n(lmax) * sizeof(CCTK_REAL));

  if (buffer == NULL)
    CCTK_ERROR("Internal error. Failed to allocate memory for the spherical "
               "harmonic array");

  return buffer;
}

CCTK_REAL Ylm(const CCTK_REAL *buffer, CCTK_INT l, CCTK_INT m,
              CCTK_REAL theta) {
  const CCTK_REAL legendre_part =
      buffer[gsl_sf_legendre_array_index(l, abs(m))];

  if (m < 0)
    return sqrt(2.0) * legendre_part * sin(-m * theta);
  else if (m > 0)
    return sqrt(2.0) * legendre_part * cos(m * theta);
  else
    return legendre_part;
}

CCTK_REAL multipolar_gaussian(CCTK_REAL *buffer, CCTK_INT lmax, CCTK_REAL x,
                              CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmx0 = x - gaussian_x0;
  const CCTK_REAL ymy0 = y - gaussian_y0;
  const CCTK_REAL zmz0 = z - gaussian_z0;

  const CCTK_REAL xmx02 = xmx0 * xmx0;
  const CCTK_REAL ymy02 = ymy0 * ymy0;
  const CCTK_REAL zmz02 = zmz0 * zmz0;

  const CCTK_REAL R2 = xmx02 + ymy02 + zmz02;
  CCTK_REAL R = sqrt(R2);

  /* TODO: Do something better for small R ?*/
  if (R < 1.0e-8)
    R = 1.0e-8;

  const CCTK_REAL theta = zmz0 / R;
  const CCTK_REAL expo = exp(-0.5 * ((R - gaussian_R0) / gaussian_sigma) *
                             ((R - gaussian_R0) / gaussian_sigma));

  /* Compute all spherical harmonics up to lmax at this point */
  const CCTK_INT ierr =
      gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, lmax, theta, -1, buffer);

  if (ierr) {
    free(buffer);
    buffer = NULL;
    CCTK_VERROR("GSL internal error: %s", gsl_strerror(ierr));
  }

  /* Compute the multipole series */
  CCTK_REAL multipole_sum =
      multipoles[0] * Ylm(buffer, 0, 0, atan2(ymy0, xmx0)) +
      multipoles[1] * Ylm(buffer, 1, -1, atan2(ymy0, xmx0)) +
      multipoles[2] * Ylm(buffer, 1, 0, atan2(ymy0, xmx0)) +
      multipoles[3] * Ylm(buffer, 1, 1, atan2(ymy0, xmx0)) +
      multipoles[4] * Ylm(buffer, 2, -2, atan2(ymy0, xmx0)) +
      multipoles[5] * Ylm(buffer, 2, -1, atan2(ymy0, xmx0)) +
      multipoles[6] * Ylm(buffer, 2, 0, atan2(ymy0, xmx0)) +
      multipoles[7] * Ylm(buffer, 2, 1, atan2(ymy0, xmx0)) +
      multipoles[8] * Ylm(buffer, 2, 2, atan2(ymy0, xmx0));

  return multipole_sum * expo;
}

void KleinGordon_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Determine which type of initial data to apply */
  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    CCTK_REAL *buffer = Ylm_buffer(2);

    CCTK_LOOP3_ALL(loop_multipolar_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      Phi[ijk] = multipolar_gaussian(buffer, 2, x[ijk], y[ijk], z[ijk]);
      K_Phi[ijk] = 0.0;
    }
    CCTK_ENDLOOP3_ALL(loop_multipolar_gaussian);

    free(buffer);

  } else if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      Phi[ijk] = exact_gaussian(0.0, x[ijk], y[ijk], z[ijk]);
      K_Phi[ijk] = dt_exact_gaussian(0.0, x[ijk], y[ijk], z[ijk]);
    }
    CCTK_ENDLOOP3_ALL(loop_exact_gaussian);
  }
}
