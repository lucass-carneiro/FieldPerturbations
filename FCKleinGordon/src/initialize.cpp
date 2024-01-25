#include "initial_conditions.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_initialize(CCTK_ARGUMENTS) {
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;
  using namespace fckg;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_initialize);
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = 0.0;

  if (CCTK_EQUALS(initial_data, "exact_gaussian")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_exact_gaussian, cctkGH, i, j, k) {

      const auto ijk{CCTK_GFINDEX3D(cctkGH, i, j, k)};

      const auto detgamma{-(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                          - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                          + gxx[ijk] * gyy[ijk] * gzz[ijk]};

      Pi[ijk] = (sqrt(detgamma)
                 * (-d_exact_gaussian_solution_dt(t, r[ijk], sigma)
                    + d_exact_gaussian_solution_dr(t, r[ijk], sigma)
                          * (betaz[ijk] * (isapprox(r[ijk], 0.0) ? 0.0 : z[ijk] / r[ijk])
                             + betay[ijk] * (isapprox(r[ijk], 0.0) ? 0.0 : y[ijk] / r[ijk])
                             + betax[ijk] * (isapprox(r[ijk], 0.0) ? 0.0 : x[ijk] / r[ijk]))))
                / alp[ijk];

      Psi_x[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma)
                   * (isapprox(r[ijk], 0.0) ? 0.0 : x[ijk] / r[ijk]);

      Psi_y[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma)
                   * (isapprox(r[ijk], 0.0) ? 0.0 : y[ijk] / r[ijk]);

      Psi_z[ijk] = d_exact_gaussian_solution_dr(t, r[ijk], sigma)
                   * (isapprox(r[ijk], 0.0) ? 0.0 : z[ijk] / r[ijk]);

      Phi[ijk] = exact_gaussian_solution(t, r[ijk], sigma);
    }
    CCTK_ENDLOOP3_ALL(loop_exact_gaussian);

  } else if (CCTK_EQUALS(initial_data, "plane_wave")) {

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_plane_wave, cctkGH, i, j, k) {
      const auto ijk{CCTK_GFINDEX3D(cctkGH, i, j, k)};

      const auto detgamma{-(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                          - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                          + gxx[ijk] * gyy[ijk] * gzz[ijk]};

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

  } else if (CCTK_EQUALS(initial_data, "gaussian")) {
    const auto amplitude{A / exact_gaussian_solution(R0, 0.0, sigma)};

#pragma omp parallel
    CCTK_LOOP3_ALL(loop_multipolar_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      const auto xL{x[ijk] - x0};
      const auto yL{y[ijk] - y0};
      const auto zL{z[ijk] - z0};
      const auto rL{sqrt(xL * xL + yL * yL + zL * zL)};

      const auto detgamma{-(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                          - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                          + gxx[ijk] * gyy[ijk] * gzz[ijk]};

      const auto sqrtg{sqrt(detgamma)};

      const auto Psi_xL{amplitude * d_exact_gaussian_solution_dr(R0, rL, sigma)
                        * (isapprox(rL, 0.0) ? 0.0 : xL / rL)};

      const auto Psi_yL{amplitude * d_exact_gaussian_solution_dr(R0, rL, sigma)
                        * (isapprox(rL, 0.0) ? 0.0 : yL / rL)};
      const auto Psi_zL{amplitude * d_exact_gaussian_solution_dr(R0, rL, sigma)
                        * (isapprox(rL, 0.0) ? 0.0 : zL / rL)};

      Phi[ijk] = amplitude * exact_gaussian_solution(R0, rL, sigma);

      Psi_x[ijk] = Psi_xL;
      Psi_y[ijk] = Psi_yL;
      Psi_z[ijk] = Psi_zL;

      Pi[ijk] = amplitude
                * (alp[ijk] / sqrtg
                   * ((betax[ijk] - (isapprox(rL, 0.0) ? 0.0 : xL / rL)) * Psi_xL
                      + (betay[ijk] - (isapprox(rL, 0.0) ? 0.0 : yL / rL)) * Psi_yL
                      + (betaz[ijk] - (isapprox(rL, 0.0) ? 0.0 : zL / rL)) * Psi_zL));
    }
    CCTK_ENDLOOP3_ALL(loop_multipolar_gaussian);
  }
}
