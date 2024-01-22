#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_zero_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_zero_rhs);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_ALL(loop_zero_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    Pi_rhs[ijk] = 0.0;
    Psi_x_rhs[ijk] = 0.0;
    Psi_y_rhs[ijk] = 0.0;
    Psi_z_rhs[ijk] = 0.0;
    Phi_rhs[ijk] = 0.0;
  }
  CCTK_ENDLOOP3_ALL(loop_zero_rhs);
}

extern "C" void FCKleinGordon_zero_flux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_zero_flux);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_ALL(loop_zero_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    F_Pi_x[ijk] = 0.0;
    F_Pi_y[ijk] = 0.0;
    F_Pi_z[ijk] = 0.0;
    F_Psi[ijk] = 0.0;
  }
  CCTK_ENDLOOP3_ALL(loop_zero_rhs);
}