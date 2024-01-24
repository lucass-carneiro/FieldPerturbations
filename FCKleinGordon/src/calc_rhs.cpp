
//clang-format off
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
//clang-format on

#include "derivatives.hpp"

#include <cmath>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_calc_rhs(CCTK_ARGUMENTS) {
  using namespace fckg;
  using std::sqrt;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_calc_rhs);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_INT(loop_rhs, cctkGH, i, j, k) {

    const auto ijk{I(cctkGH, i, j, k)};
    const deriv_data dd{i, j, k, CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)};

    const auto det_gamma{-(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
                         - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
                         + gxx[ijk] * gyy[ijk] * gzz[ijk]};

    const auto sqrtg{sqrt(det_gamma)};

    const auto S_Pi{alp[ijk] * sqrtg * field_mass * field_mass * Phi[ijk]};
    const auto S_Phi{(betax[ijk] * Psi_x[ijk] + betay[ijk] * Psi_y[ijk] + betaz[ijk] * Psi_z[ijk])
                     - alp[ijk] * Pi[ijk] / sqrtg};

    const auto dF_Pi_x_dx{global_Dx(cctkGH, dd, F_Pi_x, J11[ijk], J21[ijk], J31[ijk])};
    const auto dF_Pi_y_dy{global_Dy(cctkGH, dd, F_Pi_y, J12[ijk], J22[ijk], J32[ijk])};
    const auto dF_Pi_z_dz{global_Dz(cctkGH, dd, F_Pi_z, J13[ijk], J23[ijk], J33[ijk])};

    const auto dF_Psi_dx{global_Dx(cctkGH, dd, F_Psi, J11[ijk], J21[ijk], J31[ijk])};
    const auto dF_Psi_dy{global_Dy(cctkGH, dd, F_Psi, J12[ijk], J22[ijk], J32[ijk])};
    const auto dF_Psi_dz{global_Dz(cctkGH, dd, F_Psi, J13[ijk], J23[ijk], J33[ijk])};

    Pi_rhs[ijk] = S_Pi - (dF_Pi_x_dx + dF_Pi_y_dy + dF_Pi_z_dz);

    Psi_x_rhs[ijk] = -dF_Psi_dx;
    Psi_y_rhs[ijk] = -dF_Psi_dy;
    Psi_z_rhs[ijk] = -dF_Psi_dz;

    Phi_rhs[ijk] = S_Phi;
  }
  CCTK_ENDLOOP3_INT(loop_rhs);
}
