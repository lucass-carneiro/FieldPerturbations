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
CCTK_REAL real_spherical_harmonics(const CCTK_REAL *P_lm_array, CCTK_INT l, CCTK_INT m,
                                   CCTK_REAL phi) {

  const CCTK_REAL P_lm_cos_theta = P_lm_array[gsl_sf_legendre_array_index(l, abs(m))];

  if (m > 0)
    return P_lm_cos_theta * cos(m * phi);
  else if (m < 0)
    return P_lm_cos_theta * sin(-m * phi);
  else
    return P_lm_cos_theta;
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
CCTK_REAL multipolar_gaussian(const CCTK_REAL *multipole_array, CCTK_REAL *P_lm_array,
                              CCTK_INT lmax, CCTK_REAL R0, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                              CCTK_REAL sigma) {

  const CCTK_REAL R = sqrt(x * x + y * y + z * z);

  const CCTK_REAL cos_theta = (R < smallnes_threshold) ? z / smallnes_threshold : z / R;

  CCTK_REAL phi = 0;
  if (x > 0)
    phi = atan2(y, x);
  else if (x < 0)
    phi = atan2(y, x) + M_PI;
  else
    phi = M_PI / 2;

  // Compute all spherical harmonics up to lmax at this point
  const CCTK_INT ierr
      = gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, lmax, cos_theta, -1, P_lm_array);

  if (ierr) {
    free(P_lm_array);
    P_lm_array = NULL;
    CCTK_VERROR("GSL internal error: %s", gsl_strerror(ierr));
  }

  // Multipolar series
  const CCTK_REAL multipole_sum
      = multipole_array[0] * real_spherical_harmonics(P_lm_array, 0, 0, phi)
        + multipole_array[1] * real_spherical_harmonics(P_lm_array, 1, -1, phi)
        + multipole_array[2] * real_spherical_harmonics(P_lm_array, 1, 0, phi)
        + multipole_array[3] * real_spherical_harmonics(P_lm_array, 1, 1, phi)
        + multipole_array[4] * real_spherical_harmonics(P_lm_array, 2, -2, phi)
        + multipole_array[5] * real_spherical_harmonics(P_lm_array, 2, -1, phi)
        + multipole_array[6] * real_spherical_harmonics(P_lm_array, 2, 0, phi)
        + multipole_array[7] * real_spherical_harmonics(P_lm_array, 2, 1, phi)
        + multipole_array[8] * real_spherical_harmonics(P_lm_array, 2, 2, phi);

  return multipole_sum * base_gaussian(R - R0, sigma);
}

void KleinGordon_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ijk = 0;

  if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    const CCTK_INT max_supported_l = 2;

    CCTK_REAL *P_lm_array = create_legendre_buffer(max_supported_l);
    CCTK_REAL gaussian = 0.0, R = 0.0, gaussian_dr = 0.0, contraction = 0.0;

    CCTK_LOOP3_ALL(loop_multipolar_gaussian, cctkGH, i, j, k) {
      ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      gaussian = multipolar_gaussian(multipoles, P_lm_array, max_supported_l, gaussian_R0,
                                     x[ijk] - gaussian_x0, y[ijk] - gaussian_y0,
                                     z[ijk] - gaussian_z0, gaussian_sigma);

      R = sqrt((x[ijk] - gaussian_x0) * (x[ijk] - gaussian_x0)
               + (y[ijk] - gaussian_y0) * (y[ijk] - gaussian_y0)
               + (z[ijk] - gaussian_z0) * (z[ijk] - gaussian_z0));

      gaussian_dr = gaussian * (-(R - gaussian_R0) / (gaussian_sigma * gaussian_sigma));

      contraction = (betax[ijk] * (x[ijk] - gaussian_x0) + betay[ijk] * (y[ijk] - gaussian_y0)
                     + betaz[ijk] * (z[ijk] - gaussian_z0));

      if (contraction < 1.0e-13)
        contraction = 0.0;
      else
        contraction /= R;

      Phi[ijk] = gaussian;

      // This choice makes the gaussian move towards the origin, instead of splitting.
      K_Phi[ijk] = ((contraction - 1.0) * gaussian_dr) / (2 * alp[ijk]);
    }
    CCTK_ENDLOOP3_ALL(loop_multipolar_gaussian);

    free(P_lm_array);

  } else if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {
      ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

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

  } else if (CCTK_EQUALS(initial_data, "plane_wave")) {

    CCTK_REAL omega = 0.0;
    CCTK_REAL nx = 0.0;
    CCTK_REAL ny = 0.0;
    CCTK_REAL nz = 0.0;
    CCTK_REAL nt = 0.0;
    CCTK_REAL w = 0.0;

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_plane_wave, cctkGH, i, j, k) {
      ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      omega = sqrt(wave_number[0] * wave_number[0] + wave_number[1] * wave_number[1]
                   + wave_number[2] * wave_number[2]);
      nx = wave_number[0] * (x[ijk] - space_offset[0]);
      ny = wave_number[1] * (y[ijk] - space_offset[1]);
      nz = wave_number[2] * (z[ijk] - space_offset[2]);
      nt = -omega * time_offset;
      w = 2 * M_PI * (nx + ny + nz + nt);

      /* Since this initial data is intended to be used as a test in flat background,
       * we assume that the lapse is 1 on the equations below
       */
      Phi[ijk] = cos(w);
      K_Phi[ijk] = sin(w) * M_PI * omega;
    }
    CCTK_ENDLOOP3_ALL(loop_plane_wave);
  }
}
