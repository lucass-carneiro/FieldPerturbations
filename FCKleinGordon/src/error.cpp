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
 *  error.cpp
 *  Compare the current solution to the exact gaussian solution.
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

#include <cmath>

extern "C" void FCKleinGordon_error(CCTK_ARGUMENTS) {
  using std::fabs;
  using std::pow;
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_error);
  DECLARE_CCTK_PARAMETERS;

  /* Time values */
  const CCTK_REAL t = cctk_time;

#pragma omp parallel
  CCTK_LOOP3_INT(loop_error, cctkGH, i, j, k) {

    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    CCTK_REAL detgamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                         - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                         + gxx[ijk] * gyy[ijk] * gzz[ijk];

    const CCTK_REAL exact_Pi
        = (sqrt(detgamma)
           * (-((-1 + exp((2 * t * r[ijk]) / pow(sigma, 2))) * pow(sigma, 2)
                * (betaz[ijk] * z[ijk] / r[ijk] + y[ijk] / r[ijk] + betax[ijk] * x[ijk] / r[ijk]))
              + pow(r[ijk], 2)
                    * (-1 - exp((2 * t * r[ijk]) / pow(sigma, 2))
                       + betay[ijk] * exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2)
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

    const CCTK_REAL exact_Psi_x
        = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
            - exp((2 * t * r[ijk]) / pow(sigma, 2)) * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
           * x[ijk] / r[ijk])
          / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

    const CCTK_REAL exact_Psi_y
        = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
            - exp((2 * t * r[ijk]) / pow(sigma, 2)) * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
           * y[ijk] / r[ijk])
          / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

    const CCTK_REAL exact_Psi_z
        = ((pow(sigma, 2) + t * r[ijk] + pow(r[ijk], 2)
            - exp((2 * t * r[ijk]) / pow(sigma, 2)) * (pow(sigma, 2) - t * r[ijk] + pow(r[ijk], 2)))
           * z[ijk] / r[ijk])
          / (exp(pow(t + r[ijk], 2) / (2. * pow(sigma, 2))) * pow(sigma, 2) * pow(r[ijk], 2));

    const CCTK_REAL exact_Phi = (exp(-pow(t - r[ijk], 2) / (2. * pow(sigma, 2)))
                                 - exp(-pow(t + r[ijk], 2) / (2. * pow(sigma, 2))))
                                / r[ijk];

    Pi_error[ijk] = fabs(exact_Pi - Pi[ijk]);
    Psi_x_error[ijk] = fabs(exact_Psi_x - Psi_x[ijk]);
    Psi_y_error[ijk] = fabs(exact_Psi_y - Psi_y[ijk]);
    Psi_z_error[ijk] = fabs(exact_Psi_z - Psi_z[ijk]);
    Phi_error[ijk] = fabs(exact_Phi - Phi[ijk]);
  }
  CCTK_ENDLOOP3_INT(loop_error);
}

extern "C" void FCKleinGordon_multipatch_error(CCTK_ARGUMENTS) {
  using std::abs;
  using std::pow;
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_multipatch_error);
  DECLARE_CCTK_PARAMETERS;

  // Time values
  const CCTK_REAL t = cctk_time;

#pragma omp parallel
  CCTK_LOOP3_INT(loop_multipatch_error, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    // Minkowski backgreound is assumed for the following computations
    // clang-format off
    Pi_multipatch_error[ijk] = abs(Pi_rhs[ijk] - (-(cos(2*M_PI*((x[ijk] - space_offset[0])*wave_number[0] + (y[ijk] - space_offset[1])*wave_number[1] + (z[ijk] - space_offset[2])*wave_number[2] + (t - time_offset)*sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))*(pow(field_mass,2) + 4*pow(M_PI,2)*(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))));
    Psi_x_multipatch_error[ijk] = abs(Psi_x_rhs[ijk] - (-4*pow(M_PI,2)*cos(2*M_PI*((x[ijk] - space_offset[0])*wave_number[0] + (y[ijk] - space_offset[1])*wave_number[1] + (z[ijk] - space_offset[2])*wave_number[2] + (t - time_offset)*sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))*wave_number[0]*(-wave_number[1] + sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2)))));
    Psi_y_multipatch_error[ijk] = abs(Psi_y_rhs[ijk] - (-4*pow(M_PI,2)*cos(2*M_PI*((x[ijk] - space_offset[0])*wave_number[0] + (y[ijk] - space_offset[1])*wave_number[1] + (z[ijk] - space_offset[2])*wave_number[2] + (t - time_offset)*sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))*wave_number[1]*(-wave_number[1] + sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2)))));
    Psi_z_multipatch_error[ijk] = abs(Psi_z_rhs[ijk] - (-4*pow(M_PI,2)*cos(2*M_PI*((x[ijk] - space_offset[0])*wave_number[0] + (y[ijk] - space_offset[1])*wave_number[1] + (z[ijk] - space_offset[2])*wave_number[2] + (t - time_offset)*sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))*wave_number[2]*(-wave_number[1] + sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2)))));
    Phi_multipatch_error[ijk] = abs(Phi_rhs[ijk] - (-2*M_PI*sin(2*M_PI*((x[ijk] - space_offset[0])*wave_number[0] + (y[ijk] - space_offset[1])*wave_number[1] + (z[ijk] - space_offset[2])*wave_number[2] + (t - time_offset)*sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2))))*(-wave_number[1] + sqrt(pow(wave_number[0],2) + pow(wave_number[1],2) + pow(wave_number[2],2)))));
    // clang-format on
  }
  CCTK_ENDLOOP3_INT(loop_multipatch_error);
}