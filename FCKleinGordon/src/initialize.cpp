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

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

#include <cmath>

extern "C" void FCKleinGordon_initialize(CCTK_ARGUMENTS) {
  using std::pow;
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_initialize);
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = 0.0;

  if (CCTK_EQUALS(initial_data, "exact_gaussian")) {
#pragma omp parallel
    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {
      CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                           - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                           + gxx[ijk] * gyy[ijk] * gzz[ijk];

      if (r[ijk] < 1.0e-8) {
        Pi[ijk] = (sqrt(detgamma)
                   * (betay[ijk]
                      - (2 * (sigma - t) * (sigma + t))
                            / (exp(pow(t, 2) / (2. * pow(sigma, 2))) * pow(sigma, 4))))
                  / alp[ijk];

        Psi_x[ijk] = 0;
        Psi_y[ijk] = 0;
        Psi_z[ijk] = 0;

        Phi[ijk] = (2 * t) / (exp(pow(t, 2) / (2. * pow(sigma, 2))) * pow(sigma, 2));

      } else {
        Pi[ijk] = (sqrt(detgamma)
                   * (-((-1 + exp((2 * t * r[ijk]) / pow(sigma, 2))) * pow(sigma, 2)
                        * (betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk]
                           + betax[ijk] * x[ijk] / r[ijk]))
                      + pow(r[ijk], 2)
                            * (-1 - exp((2 * t * r[ijk]) / pow(sigma, 2))
                               + betay[ijk] * exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2)))
                                     * pow(sigma, 2)
                               + betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk]
                               + betax[ijk] * x[ijk] / r[ijk]
                               - exp((2 * t * r[ijk]) / pow(sigma, 2))
                                     * (betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk]
                                        + betax[ijk] * x[ijk] / r[ijk]))
                      + t * r[ijk]
                            * (-1 + betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk]
                               + betax[ijk] * x[ijk] / r[ijk]
                               + exp((2 * t * r[ijk]) / pow(sigma, 2))
                                     * (1 + betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk]
                                        + betax[ijk] * x[ijk] / r[ijk]))))
                  / (alp[ijk] * exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2)
                     * pow(r[ijk], 2));

        Psi_x[ijk]
            = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
                - exp((2 * t * r[ijk]) / pow(sigma, 2))
                      * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
               * x[ijk] / r[ijk])
              / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

        Psi_y[ijk]
            = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
                - exp((2 * t * r[ijk]) / pow(sigma, 2))
                      * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
               * y[ijk] / r[ijk])
              / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

        Psi_z[ijk]
            = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
                - exp((2 * t * r[ijk]) / pow(sigma, 2))
                      * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
               * z[ijk] / r[ijk])
              / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

        Phi[ijk] = (exp(-pow(t - r[ijk], 2) / (2. * pow(sigma, 2)))
                    - exp(-pow(t + r[ijk], 2) / (2. * pow(sigma, 2))))
                   / r[ijk];
      }
    }
    CCTK_ENDLOOP3_ALL(loop_exact_gaussian);

  } else if (CCTK_EQUALS(initial_data, "plane_wave")) {
#pragma omp parallel
    CCTK_LOOP3_ALL(loop_plane_wave, cctkGH, i, j, k) {
      CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                           - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                           + gxx[ijk] * gyy[ijk] * gzz[ijk];

      Pi[ijk] = (sqrt(detgamma)
                 * (betay[ijk]
                    + 2 * M_PI
                          * sin(2 * M_PI
                                * ((x[ijk] - space_offset[0]) * wave_number[0]
                                   + (y[ijk] - space_offset[1]) * wave_number[1]
                                   + (z[ijk] - space_offset[2]) * wave_number[2]
                                   + (t - time_offset)
                                         * sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                                + pow(wave_number[2], 2))))
                          * (-(betax[ijk] * wave_number[0]) - wave_number[1]
                             - betaz[ijk] * wave_number[2]
                             + sqrt(pow(wave_number[0], 2) + pow(wave_number[1], 2)
                                    + pow(wave_number[2], 2)))))
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
  }
}
