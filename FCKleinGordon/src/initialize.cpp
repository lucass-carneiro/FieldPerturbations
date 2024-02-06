#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

template <typename T> static T base_gaussian(T A, T W, T r) noexcept {
  using std::exp;
  return A * exp(-(r * r) / (W * W) / 2);
}

extern "C" void FCKleinGordon_initialize(CCTK_ARGUMENTS) {
  using std::cos;
  using std::sin;
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_initialize);
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_data, "standing_wave")) {
#pragma omp parallel
    CCTK_LOOP3_ALL(loop_stnading_wave, cctkGH, i, j, k) {

      const auto ijk{CCTK_GFINDEX3D(cctkGH, i, j, k)};

      Pi[ijk] = 0.0;

      Psi_x[ijk] = -2 * A * kx * M_PI * cos(2 * M_PI * ky * y[ijk]) * cos(2 * M_PI * kz * z[ijk])
                   * sin(2 * M_PI * kx * x[ijk]);

      Psi_y[ijk] = -2 * A * ky * M_PI * cos(2 * M_PI * kx * y[ijk]) * cos(2 * M_PI * kz * z[ijk])
                   * sin(2 * M_PI * ky * x[ijk]);

      Psi_z[ijk] = -2 * A * kz * M_PI * cos(2 * M_PI * kx * y[ijk]) * cos(2 * M_PI * ky * z[ijk])
                   * sin(2 * M_PI * kz * x[ijk]);

      Phi[ijk] = A * cos(2 * M_PI * kx * x[ijk]) * cos(2 * M_PI * ky * y[ijk])
                 * cos(2 * M_PI * kz * z[ijk]);
    }
    CCTK_ENDLOOP3_ALL(loop_stnading_wave);

  } else if (CCTK_EQUALS(initial_data, "gaussian")) {
#pragma omp parallel
    CCTK_LOOP3_ALL(loop_gaussian, cctkGH, i, j, k) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      const auto xL{x[ijk] - x0};
      const auto yL{y[ijk] - y0};
      const auto zL{z[ijk] - z0};
      const auto rL{sqrt(xL * xL + yL * yL + zL * zL)};

      const auto xL_over_rL{rL < 1.0e-3 ? 0 : xL / rL};
      const auto yL_over_rL{rL < 1.0e-3 ? 0 : yL / rL};
      const auto zL_over_rL{rL < 1.0e-3 ? 0 : zL / rL};

      const auto gaussian{base_gaussian(A, W, rL)};

      const auto detgamma{-(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                          - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                          + gxx[ijk] * gyy[ijk] * gzz[ijk]};

      const auto sqrtg{sqrt(detgamma)};

      const auto Psi_xL{xL / (W * W) * gaussian};
      const auto Psi_yL{yL / (W * W) * gaussian};
      const auto Psi_zL{zL / (W * W) * gaussian};

      Phi[ijk] = gaussian;

      Psi_x[ijk] = Psi_xL;
      Psi_y[ijk] = Psi_yL;
      Psi_z[ijk] = Psi_zL;

      Pi[ijk] = (sqrtg / alp[ijk])
                * ((betax[ijk] - xL_over_rL) * Psi_xL + (betay[ijk] - yL_over_rL) * Psi_yL
                   + (betaz[ijk] - zL_over_rL) * Psi_zL);
    }
    CCTK_ENDLOOP3_ALL(loop_gaussian);
  }
}
