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

/**
 * A constant number that represents a small value. If reals are smaller than
 * this threshold, they are taken to be zero.
 */
const CCTK_REAL smallnes_threshold = 1.0e-8;

/**
 * Computes the factorial of small integers
 *
 * @param n The number to compute the factorial.
 * @return The factorial of n
 */
CCTK_INT factorial(CCTK_INT n) {
  if (n == 0)
    return 1;
  else
    return n * factorial(n - 1);
}

/**
 * Basic Gaussian function
 *
 * @param x The argument of the gaussian
 * @param sigma The Gaussian width
 * @return The value of the Gaussian function.
 */
inline CCTK_REAL base_gaussian(CCTK_REAL x, CCTK_REAL sigma) {
  const CCTK_REAL x2 = x * x;
  const CCTK_REAL sigma2 = sigma * sigma;
  const CCTK_REAL x2_over_sigma_2 = x2 / sigma2;

  return exp(-x2_over_sigma_2 / 2);
}

/**
 * The x derivative of a basic Gaussian function
 *
 * @param x The argument of the Gaussian x derivative
 * @param sigma The gaussian width
 * @return The value of the x derivative of the gaussian function.
 */
inline CCTK_REAL base_gaussian_dx(CCTK_REAL x, CCTK_REAL sigma) {
  const CCTK_REAL sigma2 = sigma * sigma;
  return -x / sigma2 * base_gaussian(x, sigma);
}

/**
 * The general Gaussian solution of the wave equation in Minkowski
 * spacetime
 *
 * @param r The position in spherical coordiantes.
 * @param t The time which the solution is to be evaluated.
 * @param sigma The Gaussian width.
 * @return The general solution of the wave equation in Minkowski spacetime.
 */
inline CCTK_REAL gaussian_solution(CCTK_REAL r, CCTK_REAL t, CCTK_REAL sigma) {
  if (r < smallnes_threshold) {
    const CCTK_REAL sigma2 = sigma * sigma;
    return 2 * t / sigma2 * base_gaussian(t, sigma);
  } else {
    return (base_gaussian(r - t, sigma) - base_gaussian(r + t, sigma)) / r;
  }
}

/**
 * Time derivative of the general Gaussian solution of the wave equation in Minkowski
 * spacetime
 *
 * @param r The position in spherical coordiantes.
 * @param t The time which the solution derivative is to be evaluated.
 * @param sigma The Gaussian width.
 * @return The time derivative of the general solution of the wave equation in Minkowski spacetime.
 */
inline CCTK_REAL gaussian_solution_dt(CCTK_REAL r, CCTK_REAL t, CCTK_REAL sigma) {
  if (r < smallnes_threshold) {
    const CCTK_REAL sigma2 = sigma * sigma;
    const CCTK_REAL sigma4 = sigma2 * sigma2;
    const CCTK_REAL t2 = t * t;
    return 2 * (sigma2 - t2) / sigma4 * base_gaussian(t, sigma);
  } else {
    return -(base_gaussian_dx(r - t, sigma) + base_gaussian_dx(r + t, sigma)) / r;
  }
}

CCTK_REAL cartesian_gaussian_solution(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                      CCTK_REAL sigma) {
  CCTK_REAL r = sqrt(x * x + y * y + z * z);
  return gaussian_solution(r, t, sigma);
}

CCTK_REAL cartesian_gaussian_solution_dt(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                         CCTK_REAL sigma) {
  CCTK_REAL r = sqrt(x * x + y * y + z * z);
  return gaussian_solution_dt(r, t, sigma);
}

/**
 * Set all elements of a Legendre polynomials buffer array to zero
 *
 * @param buffer A pointer to the buffer array.
 * @param lmax The maximun l value that buffer contains.
 */
inline void reset_legendre_buffer(CCTK_REAL *buffer, CCTK_INT lmax) {
  for (CCTK_INT i = 0; i < gsl_sf_legendre_array_n(lmax); i++)
    buffer[i] = 0;
}

/**
 * Allocates a buffer for storing the associated legendre polynomials.
 *
 * It's assumed that the buffer ownership is moved on assignment,
 * thus the final user is responsible for freeing the memory.
 * If memory allocation fails, the function halts Cactus.
 *
 * @param lmax The maximun l value that will be used.
 * @return A pointer to a array of CCTK_REAL values that will be filled with the values of
 * the associated Legendre polynomials.
 */
CCTK_REAL *create_legendre_buffer(CCTK_INT lmax) {
  CCTK_REAL *buffer = (CCTK_REAL *)calloc(gsl_sf_legendre_array_n(lmax), sizeof(CCTK_REAL));

  if (buffer == NULL)
    CCTK_ERROR("Internal error. Failed to allocate memory for the spherical "
               "harmonic array");

  return buffer;
}

/**
 * Compute the Real spherical harmonics
 *
 * Real spherical harmonicas are defined as linear combinantions of spherical harmonics
 * or as a piecewise function on m.
 * We have implemented the last equation of
 * https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form with slight modificiations to the
 * normalization in order to reproduce the values listed here
 * https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
 *
 * @param l The l value of the spherical harmonics.
 * @param m The m value of the spherical harmonics.
 * @param P_lm_array A pointer to an array containing all values of P_{lm}(cos\theta) up
 * to the passed l value.
 * @param Phi The spherical azimuthal angle.
 * @return The value of the real spherical harmonic function.
 */
CCTK_REAL real_spherical_harmonics(CCTK_INT l, CCTK_INT m, const CCTK_REAL *P_lm_array,
                                   CCTK_REAL phi) {
  const CCTK_INT abs_m = abs(m);

  const CCTK_REAL condon_phase = pow(-1, abs_m);
  const CCTK_REAL l_minus_m_fac = factorial(l - abs_m);
  const CCTK_REAL l_plus_m_fac = factorial(l + abs_m);
  const CCTK_REAL norm
      = condon_phase * sqrt((2 * l + 1) / (4 * M_PI) * l_minus_m_fac / l_plus_m_fac);

  const CCTK_REAL P_lm_cos_theta = P_lm_array[gsl_sf_legendre_array_index(l, abs_m)];

  if (m >= 0)
    return norm * P_lm_cos_theta * cos(abs_m * phi);
  else
    return norm * P_lm_cos_theta * sin(abs_m * phi);
}

/**
 * Compute the value of the multipolar gaussian function at a given cartesian point.
 *
 * @param P_lm_array A pointer to a pre allocated buffer that will store the values of the
 * associated legendre polynomials up to lmax.
 * @param lmax The maximun value of l that will be computed.
 * @param x The x cartesian coordinate.
 * @param y The y cartesian coordiante.
 * @param z The z cartesian coordinate.
 * @return The value of the multipolar gaussian function.
 */
CCTK_REAL multipolar_gaussian(CCTK_REAL *P_lm_array, CCTK_INT lmax, CCTK_REAL x, CCTK_REAL y,
                              CCTK_REAL z, CCTK_REAL R0) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmx0 = x - gaussian_x0;
  const CCTK_REAL ymy0 = y - gaussian_y0;
  const CCTK_REAL zmz0 = z - gaussian_z0;

  const CCTK_REAL xmx02 = xmx0 * xmx0;
  const CCTK_REAL ymy02 = ymy0 * ymy0;
  const CCTK_REAL zmz02 = zmz0 * zmz0;

  const CCTK_REAL R2 = xmx02 + ymy02 + zmz02;
  CCTK_REAL R = sqrt(R2);

  const CCTK_REAL cos_theta = (R < smallnes_threshold) ? 1 : zmz0 / R;
  const CCTK_REAL phi = (R < smallnes_threshold) ? M_PI / 2 : atan2(y, x);

  // Compute all spherical harmonics up to lmax at this point
  const CCTK_INT ierr = gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, lmax, cos_theta, P_lm_array);

  if (ierr) {
    free(P_lm_array);
    P_lm_array = NULL;
    CCTK_VERROR("GSL internal error: %s", gsl_strerror(ierr));
  }

  // Multipolar series
  const CCTK_REAL multipole_sum = multipoles[0] * real_spherical_harmonics(0, 0, P_lm_array, phi)
                                  + multipoles[1] * real_spherical_harmonics(1, -1, P_lm_array, phi)
                                  + multipoles[2] * real_spherical_harmonics(1, 0, P_lm_array, phi)
                                  + multipoles[3] * real_spherical_harmonics(1, 1, P_lm_array, phi)
                                  + multipoles[4] * real_spherical_harmonics(2, -2, P_lm_array, phi)
                                  + multipoles[5] * real_spherical_harmonics(2, -1, P_lm_array, phi)
                                  + multipoles[6] * real_spherical_harmonics(2, 0, P_lm_array, phi)
                                  + multipoles[7] * real_spherical_harmonics(2, 1, P_lm_array, phi)
                                  + multipoles[8] * real_spherical_harmonics(2, 1, P_lm_array, phi);

  reset_legendre_buffer(P_lm_array, lmax);

  return multipole_sum * base_gaussian(R - R0, gaussian_sigma);
}

void KleinGordon_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // TODO: There currently is a bug in this initial condition that makes it noisy across the grid. I
  // should investigate
  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {
    const CCTK_INT max_supported_l = 2;
    CCTK_REAL *P_lm_array = create_legendre_buffer(max_supported_l);

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_multipolar_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      Phi[ijk]
          = multipolar_gaussian(P_lm_array, max_supported_l, x[ijk], y[ijk], z[ijk], gaussian_R0);
      K_Phi[ijk] = 0.0;
    }
    CCTK_ENDLOOP3_ALL(loop_multipolar_gaussian);

    free(P_lm_array);

  } else if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      /* Since this initial data represents a exact solution of the field equations
       * in a flat background, we will assume that the lapse is one and the
       * shift is zero in the equations below
       */
      Phi[ijk] = cartesian_gaussian_solution(0.0, x[ijk] - gaussian_x0, y[ijk] - gaussian_y0,
                                             z[ijk] - gaussian_z0, gaussian_sigma);
      K_Phi[ijk] = -0.5
                   * cartesian_gaussian_solution_dt(0.0, x[ijk] - gaussian_x0, y[ijk] - gaussian_y0,
                                                    z[ijk] - gaussian_z0, gaussian_sigma);
    }
    CCTK_ENDLOOP3_ALL(loop_exact_gaussian);
  }
}
