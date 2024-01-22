#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_calc_flux(CCTK_ARGUMENTS) {
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_calc_flux);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_INT(loop_rhs, cctkGH, i, j, k) {

    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    const auto det_gamma = -(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                           - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                           + gxx[ijk] * gyy[ijk] * gzz[ijk];

    const auto igxx = (-gyz[ijk] * gyz[ijk] + gyy[ijk] * gzz[ijk]) / det_gamma;
    const auto igxy = (gxz[ijk] * gyz[ijk] - gxy[ijk] * gzz[ijk]) / det_gamma;
    const auto igxz = (gxy[ijk] * gyz[ijk] - (gxz[ijk] * gyy[ijk])) / det_gamma;
    const auto igyy = (-gxz[ijk] * gxz[ijk] + gxx[ijk] * gzz[ijk]) / det_gamma;
    const auto igyz = (gxy[ijk] * gxz[ijk] - gxx[ijk] * gyz[ijk]) / det_gamma;
    const auto igzz = (-gxy[ijk] * gxy[ijk] + gxx[ijk] * gyy[ijk]) / det_gamma;

    const auto sqrtg = std::sqrt(det_gamma);

    F_Pi_x[ijk] = alp[ijk] * sqrtg * (igxx * Psi_x[ijk] + igxy * Psi_y[ijk] + igxz * Psi_z[ijk])
                  - betax[ijk] * Pi[ijk];
    F_Pi_y[ijk] = alp[ijk] * sqrtg * (igxy * Psi_x[ijk] + igyy * Psi_y[ijk] + igyz * Psi_z[ijk])
                  - betay[ijk] * Pi[ijk];
    F_Pi_z[ijk] = alp[ijk] * sqrtg * (igxz * Psi_x[ijk] + igyz * Psi_y[ijk] + igzz * Psi_z[ijk])
                  - betaz[ijk] * Pi[ijk];

    F_Psi[ijk] = alp[ijk] * Pi[ijk] / sqrtg - betax[ijk] * Psi_x[ijk] + betay[ijk] * Psi_y[ijk]
                 + betaz[ijk] * Psi_z[ijk];
  }
  CCTK_ENDLOOP3_INT(loop_rhs);
}
