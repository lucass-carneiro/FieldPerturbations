/*
 *  KleinGordonX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordonX.
 *
 *  KleinGordonX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordonX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KleinGordonX.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Initialize.cpp
 *  Fill the grid functions with initial data.
 */

#include "KleinGordonX.hpp"

/**************************
 * Std. lib. includes     *
 * and external libraries *
 **************************/
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <memory>

using namespace Loop;
using namespace std;

CCTK_REAL KleinGordonX::exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                       CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL r = sqrt((x - gaussian_x0) * (x - gaussian_x0) +
                           (y - gaussian_y0) * (y - gaussian_y0) +
                           (z - gaussian_z0) * (z - gaussian_z0));

  auto f = [&](CCTK_REAL x) { return exp(-0.5 * pow(x / gaussian_sigma, 2)); };
  auto fx = [&](CCTK_REAL x) { return -x / pow(gaussian_sigma, 2) * f(x); };

  // Use L'Hôpital's rule for small r
  if (r < 1.0e-8)
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}

CCTK_REAL KleinGordonX::dt_exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
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
 * If memory allocation fails, the function halts Cactus
 */
unique_ptr<CCTK_REAL[]> Ylm_buffer(CCTK_INT lmax) {

  auto buffer = make_unique<CCTK_REAL[]>(gsl_sf_legendre_array_n(lmax));

  if (buffer == nullptr)
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
    return pow(-1.0, -m) * legendre_part;
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
    CCTK_VWARN(CCTK_WARN_ALERT, "GSL internal error: %s", gsl_strerror(ierr));
    return 0.0;
  }

  /* Compute the multipole series */
  const CCTK_REAL multipole_sum =
      multipoles[0] * Ylm(buffer, 0, 0, atan2(ymy0, xmx0)) +
      multipoles[1] * Ylm(buffer, 1, -1, atan2(ymy0, xmx0)) +
      multipoles[2] * Ylm(buffer, 1, 0, atan2(ymy0, xmx0)) +
      multipoles[3] * Ylm(buffer, 1, 1, atan2(ymy0, xmx0)) +
      multipoles[4] * Ylm(buffer, 2, -2, atan2(ymy0, xmx0)) +
      multipoles[5] * Ylm(buffer, 2, -1, atan2(ymy0, xmx0)) +
      multipoles[6] * Ylm(buffer, 2, 0, atan2(ymy0, xmx0)) +
      multipoles[7] * Ylm(buffer, 2, 1, atan2(ymy0, xmx0)) +
      multipoles[8] * Ylm(buffer, 2, 2, atan2(ymy0, xmx0));

  const CCTK_REAL result = multipole_sum * expo;

  return result;
}

extern "C" void KleinGordonX::KleinGordonX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KleinGordonX_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_Phi(layout, Phi);
  const GF3D2<CCTK_REAL> gf_K_Phi(layout, K_Phi);

  if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

    auto exact_gaussian_lambda = [&](const PointDesc &p) {
      gf_Phi(p.I) = exact_gaussian(t, p.x, p.y, p.z);
      gf_K_Phi(p.I) = dt_exact_gaussian(t, p.x, p.y, p.z);
    };

    loop_int<0, 0, 0>(cctkGH, exact_gaussian_lambda);

  } else if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {

    auto buffer = Ylm_buffer(2);

    auto multipolar_gaussian_lambda = [&](const PointDesc &p) {
      gf_Phi(p.I) = multipolar_gaussian(buffer.get(), 2, p.x, p.y, p.z);
      gf_K_Phi(p.I) = 0.0;
    };

    loop_int<0, 0, 0>(cctkGH, multipolar_gaussian_lambda);
  }
}
