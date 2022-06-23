/*
 *  FCKleinGordon - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of FCKleinGordon.
 *
 *  FCKleinGordon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FCKleinGordon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with FCKleinGordon.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  initialize.cpp
 *  Loads the grid functions with initial data.
 */

#include "initial_conditions.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_initialize(CCTK_ARGUMENTS) {
  using std::pow;
  using std::sqrt;
  using namespace fckg;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_initialize);
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = 0.0;

  if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {

      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      const CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk])
                                 + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                                 - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                                 + gxx[ijk] * gyy[ijk] * gzz[ijk];

      Pi[ijk] = (sqrt(detgamma)
                 * (-d_exact_gaussian_solution_dt(t, r[ijk], sigma)
                    + d_exact_gaussian_solution_dr(t, r[ijk], sigma)
                          * (betaz[ijk] * z[ijk] / r[ijk] + betay[ijk] * y[ijk] / r[ijk]
                             + betax[ijk] * x[ijk] / r[ijk])))
                / alp[ijk];

      Psi_x[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma) * x[ijk] / r[ijk];

      Psi_y[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma) * y[ijk] / r[ijk];

      Psi_z[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma) * z[ijk] / r[ijk];

      Phi[ijk] = exact_gaussian_solution(t, r[ijk], sigma);
    }
    CCTK_ENDLOOP3_ALL(loop_exact_gaussian);

  } else if (CCTK_EQUALS(initial_data, "plane_wave")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_plane_wave, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      const CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk])
                                 + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                                 - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                                 + gxx[ijk] * gyy[ijk] * gzz[ijk];

      Pi[ijk]
          = (2 * sqrt(detgamma) * M_PI
             * sin(2 * M_PI
                   * ((x[ijk] - space_offset[0]) * wave_number[0]
                      + (y[ijk] - space_offset[1]) * wave_number[1]
                      + (z[ijk] - space_offset[2]) * wave_number[2]
                      + (t - time_offset)
                            * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                   + pow(wave_number[2], 2))))
             * (-(betax[ijk] * wave_number[0]) - betay[ijk] * wave_number[1]
                - betaz[ijk] * wave_number[2]
                + sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2) + pow(wave_number[2], 2))))
            / alp[ijk];

      Psi_x[ijk] = -2 * M_PI
                   * sin(2 * M_PI
                         * ((x[ijk] - space_offset[0]) * wave_number[0]
                            + (y[ijk] - space_offset[1]) * wave_number[1]
                            + (z[ijk] - space_offset[2]) * wave_number[2]
                            + (t - time_offset)
                                  * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                         + pow(wave_number[2], 2))))
                   * wave_number[0];

      Psi_y[ijk] = -2 * M_PI
                   * sin(2 * M_PI
                         * ((x[ijk] - space_offset[0]) * wave_number[0]
                            + (y[ijk] - space_offset[1]) * wave_number[1]
                            + (z[ijk] - space_offset[2]) * wave_number[2]
                            + (t - time_offset)
                                  * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                         + pow(wave_number[2], 2))))
                   * wave_number[1];

      Psi_z[ijk] = -2 * M_PI
                   * sin(2 * M_PI
                         * ((x[ijk] - space_offset[0]) * wave_number[0]
                            + (y[ijk] - space_offset[1]) * wave_number[1]
                            + (z[ijk] - space_offset[2]) * wave_number[2]
                            + (t - time_offset)
                                  * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                         + pow(wave_number[2], 2))))
                   * wave_number[2];

      Phi[ijk] = cos(2 * M_PI
                     * ((x[ijk] - space_offset[0]) * wave_number[0]
                        + (y[ijk] - space_offset[1]) * wave_number[1]
                        + (z[ijk] - space_offset[2]) * wave_number[2]
                        + (t - time_offset)
                              * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                     + pow(wave_number[2], 2))));
    }
    CCTK_ENDLOOP3_ALL(loop_plane_wave);

  } else if (CCTK_EQUALS(initial_data, "multipolar_gaussian")) {
#pragma omp parallel
    CCTK_LOOP3_ALL(loop_multipolar_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      const CCTK_REAL xL = x[ijk] - x0;
      const CCTK_REAL yL = y[ijk] - y0;
      const CCTK_REAL zL = z[ijk] - z0;
      const CCTK_REAL rL = sqrt(xL * xL + yL * yL + zL * zL);

      const CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk])
                                 + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                                 - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                                 + gxx[ijk] * gyy[ijk] * gzz[ijk];

      Phi[ijk] = multipole_sum(multipoles, xL, yL, zL) * base_gaussian(rL - R0, sigma);

      const CCTK_REAL Psi_xL
          = d_multipole_sum_dx(multipoles, xL, yL, zL) * base_gaussian(rL - R0, sigma)
            + multipole_sum(multipoles, xL, yL, zL) * d_base_gaussian_dr(rL - R0, sigma)
                  * (isapprox(rL, 0.0) ? 0.0 : xL / rL);

      const CCTK_REAL Psi_yL
          = d_multipole_sum_dy(multipoles, xL, yL, zL) * base_gaussian(rL - R0, sigma)
            + multipole_sum(multipoles, xL, yL, zL) * d_base_gaussian_dr(rL - R0, sigma)
                  * (isapprox(rL, 0.0) ? 0.0 : yL / rL);

      const CCTK_REAL Psi_zL
          = d_multipole_sum_dz(multipoles, xL, yL, zL) * base_gaussian(rL - R0, sigma)
            + multipole_sum(multipoles, xL, yL, zL) * d_base_gaussian_dr(rL - R0, sigma)
                  * (isapprox(rL, 0.0) ? 0.0 : zL / rL);

      Psi_x[ijk] = Psi_xL;
      Psi_y[ijk] = Psi_yL;
      Psi_z[ijk] = Psi_zL;

      // This choice of Pi makes the gaussian be infalling
      Pi[ijk] = sqrt(detgamma) / alp[ijk]
                * ((betax[ijk] - (isapprox(rL, 0.0) ? 0.0 : xL / rL)) * Psi_xL
                   + (betay[ijk] - (isapprox(rL, 0.0) ? 0.0 : yL / rL)) * Psi_yL
                   + (betaz[ijk] - (isapprox(rL, 0.0) ? 0.0 : zL / rL)) * Psi_zL);
    }
    CCTK_ENDLOOP3_ALL(loop_multipolar_gaussian);
  }
}
