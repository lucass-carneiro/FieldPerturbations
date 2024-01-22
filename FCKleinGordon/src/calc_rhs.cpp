#include "derivatives.hpp"

//clang-format off
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
//clang-format on

#include <cmath>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_calc_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_calc_rhs);
  DECLARE_CCTK_PARAMETERS;
  DECLARE_DERIVATIVE_FACTORS_4;

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

    const auto S_Pi = alp[ijk] * sqrtg * field_mass * field_mass * Phi[ijk];
    const auto S_Phi = betax[ijk] * Psi_x[ijk] + betay[ijk] * Psi_y[ijk] + betaz[ijk] * Psi_z[ijk]
                       - alp[ijk] * Pi[ijk] / sqrtg;

    const auto dF_Pi_x_dx = global_Dx(4, F_Pi_x);
    const auto dF_Pi_y_dy = global_Dy(4, F_Pi_y);
    const auto dF_Pi_z_dz = global_Dz(4, F_Pi_z);

    const auto dF_Psi_dx = global_Dz(4, F_Psi);
    const auto dF_Psi_dy = global_Dz(4, F_Psi);
    const auto dF_Psi_dz = global_Dz(4, F_Psi);

    Pi_rhs[ijk] = S_Pi - (dF_Pi_x_dx + dF_Pi_y_dy + dF_Pi_z_dz);

    Psi_x_rhs[ijk] = -dF_Psi_dx;
    Psi_y_rhs[ijk] = -dF_Psi_dy;
    Psi_z_rhs[ijk] = -dF_Psi_dz;

    Phi_rhs[ijk] = S_Phi;
  }
  CCTK_ENDLOOP3_INT(loop_rhs);
}
